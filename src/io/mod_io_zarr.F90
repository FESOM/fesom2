module mod_io_zarr
    ! Hand-rolled Zarr v2 writer for FESOM3 model output (M9, Task 0.1/0.2).
    !
    ! A Zarr v2 store is "just files" (spec: https://zarr.readthedocs.io/en/stable/spec/v2.html):
    !   store.zarr/
    !     .zgroup           {"zarr_format": 2}
    !     .zattrs           { ...group attrs... }
    !     <var>/
    !       .zarray         {"zarr_format":2,"shape":[..],"chunks":[..],"dtype":"<f4|<f8|<i4",
    !                        "compressor":null|{"id":"lz4",..},"fill_value":<num|null>,
    !                        "order":"C","filters":null,"dimension_separator":"."}
    !       .zattrs         {"_ARRAY_DIMENSIONS":["time","nod2"], "units":.., "long_name":.., ..}
    !       0, 0.1, 1.0 ..  chunk files (raw little-endian bytes; codec none here, lz4 in Task 0.3)
    !
    ! Design notes (cited in the M9 plan):
    !  - C-ORDER gotcha: Zarr stores arrays C row-major; Fortran is column-major, so a chunk is
    !    written transposed (the LAST Zarr dim varies fastest on disk). The caller passes a Fortran
    !    array whose i-th index is Zarr dim i (dim1 = shape(1) = slowest); we pack to C order.
    !  - LAST/PARTIAL chunk: every chunk file is full chunk-size on disk; the final global chunk
    !    along a dim is padded to full size with fill_value (the reader truncates to `shape`). This
    !    is the ONLY padding (the canonical write scheme has no per-rank padding).
    !  - Little-endian: dtype strings are "<.." and we write native bytes; x86-64 is LE (asserted by
    !    the round-trip gate via Python which honours the "<" byte-order).
    !
    ! Style mirrors src/io/mod_io_netcdf.F90 (zarr_check(ok,ctx) -> error stop) and the project
    ! conventions (WP data, explicit public list, lazy/typed handles).
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int8, error_unit
    use, intrinsic :: iso_c_binding,   only: c_ptr, c_loc, c_f_pointer, c_int, c_char, c_null_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    use mod_precision, only: WP
    implicit none
    private

    interface
        ! int mkdir(const char *pathname, mode_t mode);  — POSIX, used instead of
        ! execute_command_line so we never fork (fork after MPI_Init segfaults on Levante's
        ! OpenMPI/vader stack; the project's MPI-fork gotcha). Returns 0 / -1 (EEXIST ignored).
        function c_mkdir(path, mode) bind(C, name="mkdir") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path
            integer(c_int), value :: mode
            integer(c_int)        :: r
        end function c_mkdir
    end interface

#ifdef HAVE_LZ4
    interface
        ! int LZ4_compressBound(int inputSize);
        function c_lz4_bound(insz) bind(C, name="LZ4_compressBound") result(r)
            import :: c_int
            integer(c_int), value :: insz
            integer(c_int)        :: r
        end function c_lz4_bound
        ! int LZ4_compress_default(const char* src, char* dst, int srcSize, int dstCapacity);
        function c_lz4_compress(src, dst, srcsz, dstcap) bind(C, name="LZ4_compress_default") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value    :: src, dst
            integer(c_int), value :: srcsz, dstcap
            integer(c_int)        :: r
        end function c_lz4_compress
        ! int LZ4_decompress_safe(const char* src, char* dst, int compressedSize, int dstCapacity);
        ! the inverse of LZ4_compress_default — decodes a chunk for the chunk_time>1 read-modify-write.
        function c_lz4_decompress(src, dst, csize, dstcap) bind(C, name="LZ4_decompress_safe") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value    :: src, dst
            integer(c_int), value :: csize, dstcap
            integer(c_int)        :: r
        end function c_lz4_decompress
    end interface
#endif

    integer, parameter, public :: ZARR_MAXDIM = 4

    ! A store = a root directory holding a Zarr group. `meta` accumulates the consolidated-metadata
    ! members ("key":<json-object>) as arrays are defined, so zarr_consolidate can emit .zmetadata.
    type, public :: t_zarr_store
        character(len=4096)           :: path  = ''
        character(len=:), allocatable :: meta
        integer                       :: nmeta = 0
    end type t_zarr_store

    ! One array (= one variable) within a store.
    type, public :: t_zarr_array
        character(len=64) :: name  = ''
        integer           :: ndim  = 0
        integer           :: dims(ZARR_MAXDIM)   = 0   ! Zarr/C dim order (dims(1) slowest)
        integer           :: chunks(ZARR_MAXDIM) = 0
        character(len=8)  :: dtype = '<f8'             ! '<f8' | '<f4' | '<i4'
        character(len=8)  :: codec = 'none'            ! 'none' | 'lz4' (Task 0.3)
        logical           :: has_fill = .true.
        real(real64)      :: fill  = 0.0_real64        ! formatted per dtype on write
    end type t_zarr_array

    ! A growing JSON object body (comma-joined "key":value members), built by the zattr_* helpers
    ! and written by zarr_create_store / zarr_define_array.
    type, public :: t_zarr_attrs
        character(len=:), allocatable :: body
        integer :: n = 0
    end type t_zarr_attrs

    public :: zarr_create_store, zarr_define_array, zarr_consolidate, zarr_rewrite_zarray
    public :: zarr_array_init, zarr_attrs_init
    public :: zattr_str, zattr_int, zattr_real
    public :: zattr_str_arr, zattr_int_arr, zattr_real_arr
    public :: zarr_write_whole, zarr_write_chunk, zarr_read_chunk
    public :: zarr_check, zarr_mkdir

    interface zarr_write_whole
        module procedure zarr_write_whole_1d_real, zarr_write_whole_2d_real, &
                         zarr_write_whole_1d_int,  zarr_write_whole_2d_int
    end interface zarr_write_whole

    interface zarr_write_chunk
        module procedure zarr_write_chunk_1d_real, zarr_write_chunk_2d_real, &
                         zarr_write_chunk_3d_real, &
                         zarr_write_chunk_1d_int,  zarr_write_chunk_2d_int
    end interface zarr_write_chunk

    ! Read a chunk file back into a FULL chunk-sized Fortran array (codec-decoded, C->Fortran
    ! un-transpose). A missing file => the array filled with arr%fill. Used by mod_io_means for the
    ! chunk_time>1 read-modify-write append (Task 2.6): the partial time-chunk already on disk is
    ! read, the new time slot is overwritten, and the whole chunk is re-written.
    interface zarr_read_chunk
        module procedure zarr_read_chunk_1d_real, zarr_read_chunk_2d_real, zarr_read_chunk_3d_real
    end interface zarr_read_chunk

contains

    ! ------------------------------------------------------------------ checks / io primitives

    subroutine zarr_check(ok, ctx)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: ctx
        if (.not. ok) then
            write(error_unit,*) '[mod_io_zarr] error: ', trim(ctx)
            error stop 1
        end if
    end subroutine zarr_check

    ! Recursive `mkdir -p` via POSIX mkdir (no fork — safe after MPI_Init). Creates each path
    ! component; EEXIST (-1) is ignored. A truly-uncreatable dir surfaces later as a failed chunk
    ! open (error stop). Rank-0-only + pre-barrier at multi-rank (the store-create ordering).
    subroutine zarr_mkdir(path)
        character(len=*), intent(in) :: path
        integer :: i, n, r
        n = len_trim(path)
        do i = 2, n
            if (path(i:i) == '/') r = c_mkdir(path(1:i-1)//c_null_char, int(493, c_int))  ! 0755
        end do
        if (n >= 1) r = c_mkdir(path(1:n)//c_null_char, int(493, c_int))
    end subroutine zarr_mkdir

    ! Write a text blob (e.g. a JSON document) verbatim + a trailing newline, via stream so any
    ! length is safe (no formatted-record truncation).
    subroutine write_text_file(path, text)
        character(len=*), intent(in) :: path
        character(len=*), intent(in) :: text
        integer :: u, ios
        open(newunit=u, file=trim(path), status='replace', action='write', &
             form='unformatted', access='stream', iostat=ios)
        call zarr_check(ios == 0, 'open(write) '//trim(path))
        write(u) text//char(10)
        close(u)
    end subroutine write_text_file

    ! ------------------------------------------------------------------ JSON value formatting

    ! Minimal JSON string escaping (\ and ").
    function json_escape(s) result(o)
        character(len=*), intent(in)  :: s
        character(len=:), allocatable :: o
        integer :: i
        o = ''
        do i = 1, len_trim(s)
            select case (s(i:i))
            case ('\');  o = o // '\\'
            case ('"');  o = o // '\"'
            case default; o = o // s(i:i)
            end select
        end do
    end function json_escape

    function json_quote(s) result(o)
        character(len=*), intent(in)  :: s
        character(len=:), allocatable :: o
        o = '"'//json_escape(s)//'"'
    end function json_quote

    function json_int(i) result(o)
        integer, intent(in)           :: i
        character(len=:), allocatable :: o
        character(len=32) :: tmp
        write(tmp, '(i0)') i
        o = trim(tmp)
    end function json_int

    ! real64 -> JSON number string with round-trip precision. NaN/Inf use Python's json-accepted
    ! tokens (json.loads parses NaN/Infinity/-Infinity by default; zarr writes these for float fill).
    function json_real(r) result(o)
        real(real64), intent(in)      :: r
        character(len=:), allocatable :: o
        character(len=32) :: tmp
        if (ieee_is_nan(r)) then
            o = 'NaN'
        else if (.not. ieee_is_finite(r)) then
            if (r > 0.0_real64) then; o = 'Infinity'; else; o = '-Infinity'; end if
        else
            write(tmp, '(es24.16)') r          ! 17 sig figs, e.g. 1.2340000000000000E+02
            o = trim(adjustl(tmp))
        end if
    end function json_real

    ! ------------------------------------------------------------------ attrs builder

    subroutine zarr_attrs_init(a)
        type(t_zarr_attrs), intent(inout) :: a
        a%body = ''
        a%n    = 0
    end subroutine zarr_attrs_init

    subroutine attr_append(a, member)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: member
        if (.not. allocated(a%body)) a%body = ''
        if (a%n > 0) a%body = a%body // ','
        a%body = a%body // member
        a%n = a%n + 1
    end subroutine attr_append

    subroutine zattr_str(a, key, val)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key, val
        call attr_append(a, json_quote(key)//':'//json_quote(val))
    end subroutine zattr_str

    subroutine zattr_int(a, key, val)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key
        integer,            intent(in)    :: val
        call attr_append(a, json_quote(key)//':'//json_int(val))
    end subroutine zattr_int

    subroutine zattr_real(a, key, val)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key
        real(real64),       intent(in)    :: val
        call attr_append(a, json_quote(key)//':'//json_real(val))
    end subroutine zattr_real

    ! 1-D array of strings, e.g. _ARRAY_DIMENSIONS. Pass same-length elements
    ! ([character(len=8) :: 'time','nz1','nod2']); each is trimmed.
    subroutine zattr_str_arr(a, key, vals)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key
        character(len=*),   intent(in)    :: vals(:)
        character(len=:), allocatable :: lst
        integer :: i
        lst = '['
        do i = 1, size(vals)
            if (i > 1) lst = lst // ','
            lst = lst // json_quote(trim(vals(i)))
        end do
        lst = lst // ']'
        call attr_append(a, json_quote(key)//':'//lst)
    end subroutine zattr_str_arr

    subroutine zattr_int_arr(a, key, vals)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key
        integer,            intent(in)    :: vals(:)
        character(len=:), allocatable :: lst
        integer :: i
        lst = '['
        do i = 1, size(vals)
            if (i > 1) lst = lst // ','
            lst = lst // json_int(vals(i))
        end do
        lst = lst // ']'
        call attr_append(a, json_quote(key)//':'//lst)
    end subroutine zattr_int_arr

    subroutine zattr_real_arr(a, key, vals)
        type(t_zarr_attrs), intent(inout) :: a
        character(len=*),   intent(in)    :: key
        real(real64),       intent(in)    :: vals(:)
        character(len=:), allocatable :: lst
        integer :: i
        lst = '['
        do i = 1, size(vals)
            if (i > 1) lst = lst // ','
            lst = lst // json_real(vals(i))
        end do
        lst = lst // ']'
        call attr_append(a, json_quote(key)//':'//lst)
    end subroutine zattr_real_arr

    ! ------------------------------------------------------------------ store / array definition

    ! Append one consolidated-metadata member ("key": <json-object>) to the store accumulator.
    subroutine meta_append(store, key, jsonval)
        type(t_zarr_store), intent(inout) :: store
        character(len=*),   intent(in)    :: key, jsonval
        if (.not. allocated(store%meta)) store%meta = ''
        if (store%nmeta > 0) store%meta = store%meta // ','
        store%meta = store%meta // json_quote(key) // ':' // jsonval
        store%nmeta = store%nmeta + 1
    end subroutine meta_append

    ! Create the store root dir + .zgroup, and the root .zattrs (optional attrs).
    subroutine zarr_create_store(store, path, attrs)
        type(t_zarr_store), intent(out)          :: store
        character(len=*),   intent(in)           :: path
        type(t_zarr_attrs), intent(in), optional :: attrs
        character(len=:), allocatable :: body
        store%path  = path
        store%meta  = ''
        store%nmeta = 0
        call zarr_mkdir(trim(store%path))
        call write_text_file(trim(store%path)//'/.zgroup', '{"zarr_format": 2}')
        body = '{}'
        if (present(attrs)) then
            if (allocated(attrs%body)) body = '{'//attrs%body//'}'
        end if
        call write_text_file(trim(store%path)//'/.zattrs', body)
        call meta_append(store, '.zgroup', '{"zarr_format": 2}')
        call meta_append(store, '.zattrs', body)
    end subroutine zarr_create_store

    ! Write the optional consolidated metadata (.zmetadata) — concatenation of every .zgroup/
    ! .zarray/.zattrs as a single JSON doc, for fast opens (xr.open_zarr(consolidated=True)).
    ! Call on the metadata-writing rank (rank 0) after all arrays are defined.
    subroutine zarr_consolidate(store)
        type(t_zarr_store), intent(in) :: store
        character(len=:), allocatable :: m
        m = ''
        if (allocated(store%meta)) m = store%meta
        call write_text_file(trim(store%path)//'/.zmetadata', &
             '{"zarr_consolidated_format":1,"metadata":{'//m//'}}')
    end subroutine zarr_consolidate

    subroutine zarr_array_init(arr, name, dims, chunks, dtype, fill, has_fill, codec)
        type(t_zarr_array), intent(out)          :: arr
        character(len=*),   intent(in)           :: name
        integer,            intent(in)           :: dims(:), chunks(:)
        character(len=*),   intent(in)           :: dtype
        real(real64),       intent(in), optional :: fill
        logical,            intent(in), optional :: has_fill
        character(len=*),   intent(in), optional :: codec
        call zarr_check(size(dims) == size(chunks), 'array_init dims/chunks rank mismatch '//trim(name))
        call zarr_check(size(dims) <= ZARR_MAXDIM,  'array_init ndim > ZARR_MAXDIM '//trim(name))
        arr%name  = name
        arr%ndim  = size(dims)
        arr%dims  = 0
        arr%chunks= 0
        arr%dims(1:arr%ndim)   = dims
        arr%chunks(1:arr%ndim) = chunks
        arr%dtype = dtype
        arr%codec = 'none'
        if (present(codec))    arr%codec    = codec
        arr%has_fill = .true.
        if (present(has_fill)) arr%has_fill = has_fill
        arr%fill  = 0.0_real64
        if (present(fill))     arr%fill     = fill
    end subroutine zarr_array_init

    function json_fill(arr) result(o)
        type(t_zarr_array), intent(in) :: arr
        character(len=:), allocatable  :: o
        if (.not. arr%has_fill) then
            o = 'null'
        else if (arr%dtype(1:2) == '<i') then
            o = json_int(nint(arr%fill))
        else
            o = json_real(arr%fill)
        end if
    end function json_fill

    function json_compressor(arr) result(o)
        type(t_zarr_array), intent(in) :: arr
        character(len=:), allocatable  :: o
        select case (trim(arr%codec))
        case ('none'); o = 'null'
        case ('lz4');  o = '{"id":"lz4","acceleration":1}'
        case default;  call zarr_check(.false., 'unknown codec '//trim(arr%codec)); o = 'null'
        end select
    end function json_compressor

    ! Serialize an array's .zarray document (full v2 schema). Reads arr%dims/chunks live, so a
    ! growing array (time dim) re-serializes with its current shape via zarr_rewrite_zarray.
    function zarray_json(arr) result(zarray)
        type(t_zarr_array), intent(in) :: arr
        character(len=:), allocatable  :: zarray
        integer :: d
        zarray = '{"zarr_format":2,"shape":['
        do d = 1, arr%ndim
            if (d > 1) zarray = zarray // ','
            zarray = zarray // json_int(arr%dims(d))
        end do
        zarray = zarray // '],"chunks":['
        do d = 1, arr%ndim
            if (d > 1) zarray = zarray // ','
            zarray = zarray // json_int(arr%chunks(d))
        end do
        zarray = zarray // '],"dtype":'//json_quote(trim(arr%dtype))// &
                 ',"compressor":'//json_compressor(arr)// &
                 ',"fill_value":'//json_fill(arr)// &
                 ',"order":"C","filters":null,"dimension_separator":"."}'
    end function zarray_json

    ! Make the array subdir + write .zarray (full v2 schema) and .zattrs.
    subroutine zarr_define_array(store, arr, attrs)
        type(t_zarr_store), intent(inout)        :: store
        type(t_zarr_array), intent(in)           :: arr
        type(t_zarr_attrs), intent(in), optional :: attrs
        character(len=:), allocatable :: dir, zarray, body
        dir = trim(store%path)//'/'//trim(arr%name)
        call zarr_mkdir(dir)
        zarray = zarray_json(arr)
        call write_text_file(dir//'/.zarray', zarray)
        body = '{}'
        if (present(attrs)) then
            if (allocated(attrs%body)) body = '{'//attrs%body//'}'
        end if
        call write_text_file(dir//'/.zattrs', body)
        call meta_append(store, trim(arr%name)//'/.zarray', zarray)
        call meta_append(store, trim(arr%name)//'/.zattrs', body)
    end subroutine zarr_define_array

    ! Re-write ONLY <var>/.zarray from the (possibly mutated) handle — for a growing array whose
    ! time dim shape[0] is bumped each output record (chunk_time=1 append: new chunk files + this
    ! shape bump, NO read-modify-write). No mkdir / .zattrs / consolidation-meta touch. Rank 0 only.
    subroutine zarr_rewrite_zarray(store, arr)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        call write_text_file(trim(store%path)//'/'//trim(arr%name)//'/.zarray', zarray_json(arr))
    end subroutine zarr_rewrite_zarray

    ! ------------------------------------------------------------------ chunk path + raw writers

    function chunk_path(store, arr, cidx) result(p)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        character(len=:), allocatable  :: p, nm
        integer :: d
        nm = ''
        do d = 1, size(cidx)
            if (d > 1) nm = nm // '.'
            nm = nm // json_int(cidx(d))
        end do
        p = trim(store%path)//'/'//trim(arr%name)//'/'//nm
    end function chunk_path

    ! Write one chunk's raw little-endian byte stream, applying the codec.
    !   none : raw bytes verbatim.
    !   lz4  : numcodecs framing = 4-byte LE int32 decompressed length + LZ4 block
    !          (LZ4_compress_default == LZ4_compress_fast acceleration=1, matching .zarray).
    subroutine write_chunk_bytes(path, raw, codec)
        character(len=*), intent(in)                  :: path
        integer(int8),    intent(in), target, contiguous :: raw(:)
        character(len=*), intent(in)                  :: codec
        integer :: u, ios
        open(newunit=u, file=trim(path), status='replace', action='write', &
             form='unformatted', access='stream', iostat=ios)
        call zarr_check(ios == 0, 'open(chunk) '//trim(path))
        select case (trim(codec))
        case ('none')
            write(u) raw
        case ('lz4')
#ifdef HAVE_LZ4
            block
                integer(int8), allocatable, target :: comp(:)
                integer(c_int) :: nbytes, bound, csize
                nbytes = int(size(raw), c_int)
                bound  = c_lz4_bound(nbytes)
                allocate(comp(bound))
                csize  = c_lz4_compress(c_loc(raw(1)), c_loc(comp(1)), nbytes, bound)
                call zarr_check(csize > 0, 'LZ4_compress_default failed '//trim(path))
                write(u) int(size(raw), int32)        ! 4-byte LE decompressed length
                write(u) comp(1:csize)
            end block
#else
            call zarr_check(.false., 'lz4 codec requested but HAVE_LZ4 off (relink with liblz4)')
#endif
        case default
            call zarr_check(.false., 'unknown codec '//trim(codec))
        end select
        close(u)
    end subroutine write_chunk_bytes

    ! Typed buffer -> raw little-endian bytes (bit-copy; x86-64 native = LE) -> chunk file.
    ! Byte widths are fixed by the explicit kinds (real64=8, real32=4, int32=4), independent of WP.
    ! The bytes are a pointer VIEW of buf (c_f_pointer), not a copy: transfer() built two
    ! chunk-size temporaries that ifort places on the STACK (2 x 184 MB for a 46 x 500000 <f8
    ! restart chunk of ORCA025 -> SIGSEGV under the job's ulimit -s 204800).
    subroutine emit_r8(path, buf, codec)
        character(len=*), intent(in) :: path, codec
        real(real64),     intent(in), target, contiguous :: buf(:)
        integer(int8), pointer, contiguous :: raw(:)
        call c_f_pointer(c_loc(buf(1)), raw, [8 * size(buf)])
        call write_chunk_bytes(path, raw, codec)
    end subroutine emit_r8

    subroutine emit_r4(path, buf, codec)
        character(len=*), intent(in) :: path, codec
        real(real32),     intent(in), target, contiguous :: buf(:)
        integer(int8), pointer, contiguous :: raw(:)
        call c_f_pointer(c_loc(buf(1)), raw, [4 * size(buf)])
        call write_chunk_bytes(path, raw, codec)
    end subroutine emit_r4

    subroutine emit_i4(path, buf, codec)
        character(len=*), intent(in) :: path, codec
        integer(int32),   intent(in), target, contiguous :: buf(:)
        integer(int8), pointer, contiguous :: raw(:)
        call c_f_pointer(c_loc(buf(1)), raw, [4 * size(buf)])
        call write_chunk_bytes(path, raw, codec)
    end subroutine emit_i4

    ! ------------------------------------------------------------------ chunk writers (real input)

    ! 1-D real chunk. data(1:n0) maps to Zarr indices [cidx(1)*c0 .. +n0-1]; padded to c0 with fill.
    subroutine zarr_write_chunk_1d_real(store, arr, cidx, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        real(real64),       intent(in) :: data(:)
        real(real64), allocatable :: b8(:)
        real(real32), allocatable :: b4(:)
        integer :: c0, n0, i
        c0 = arr%chunks(1); n0 = size(data)
        if (trim(arr%dtype) == '<f8') then
            allocate(b8(c0)); b8 = arr%fill
            do i = 1, n0; b8(i) = data(i); end do
            call emit_r8(chunk_path(store, arr, cidx), b8, trim(arr%codec))
        else if (trim(arr%dtype) == '<f4') then
            allocate(b4(c0)); b4 = real(arr%fill, real32)
            do i = 1, n0; b4(i) = real(data(i), real32); end do
            call emit_r4(chunk_path(store, arr, cidx), b4, trim(arr%codec))
        else
            call zarr_check(.false., 'real data into non-float array '//trim(arr%name))
        end if
    end subroutine zarr_write_chunk_1d_real

    ! 2-D real chunk. data(i0,i1): i0=Zarr dim1 (slow), i1=Zarr dim2 (fast). Packed to C order
    ! (dim2 fastest) and padded to chunks(1)xchunks(2) with fill.
    subroutine zarr_write_chunk_2d_real(store, arr, cidx, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        real(real64),       intent(in) :: data(:,:)
        real(real64), allocatable :: b8(:)
        real(real32), allocatable :: b4(:)
        integer :: c0, c1, n0, n1, i0, i1
        c0 = arr%chunks(1); c1 = arr%chunks(2)
        n0 = size(data,1);  n1 = size(data,2)
        if (trim(arr%dtype) == '<f8') then
            allocate(b8(c0*c1)); b8 = arr%fill
            do i0 = 1, n0
                do i1 = 1, n1
                    b8((i0-1)*c1 + i1) = data(i0,i1)
                end do
            end do
            call emit_r8(chunk_path(store, arr, cidx), b8, trim(arr%codec))
        else if (trim(arr%dtype) == '<f4') then
            allocate(b4(c0*c1)); b4 = real(arr%fill, real32)
            do i0 = 1, n0
                do i1 = 1, n1
                    b4((i0-1)*c1 + i1) = real(data(i0,i1), real32)
                end do
            end do
            call emit_r4(chunk_path(store, arr, cidx), b4, trim(arr%codec))
        else
            call zarr_check(.false., 'real data into non-float array '//trim(arr%name))
        end if
    end subroutine zarr_write_chunk_2d_real

    ! 3-D real chunk. data(i0,i1,i2): i0=Zarr dim1 (slow), i1=dim2, i2=dim3 (fast). Packed to C order
    ! (dim3 fastest) and padded to chunks(1)xchunks(2)xchunks(3) with fill. Used for field output
    ! (time, nz, nod2): the per-record chunk is (1, nlev, C).
    subroutine zarr_write_chunk_3d_real(store, arr, cidx, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        real(real64),       intent(in) :: data(:,:,:)
        real(real64), allocatable :: b8(:)
        real(real32), allocatable :: b4(:)
        integer :: c0, c1, c2, n0, n1, n2, i0, i1, i2
        c0 = arr%chunks(1); c1 = arr%chunks(2); c2 = arr%chunks(3)
        n0 = size(data,1);  n1 = size(data,2);  n2 = size(data,3)
        if (trim(arr%dtype) == '<f8') then
            allocate(b8(c0*c1*c2)); b8 = arr%fill
            do i0 = 1, n0
                do i1 = 1, n1
                    do i2 = 1, n2
                        b8((i0-1)*c1*c2 + (i1-1)*c2 + i2) = data(i0,i1,i2)
                    end do
                end do
            end do
            call emit_r8(chunk_path(store, arr, cidx), b8, trim(arr%codec))
        else if (trim(arr%dtype) == '<f4') then
            allocate(b4(c0*c1*c2)); b4 = real(arr%fill, real32)
            do i0 = 1, n0
                do i1 = 1, n1
                    do i2 = 1, n2
                        b4((i0-1)*c1*c2 + (i1-1)*c2 + i2) = real(data(i0,i1,i2), real32)
                    end do
                end do
            end do
            call emit_r4(chunk_path(store, arr, cidx), b4, trim(arr%codec))
        else
            call zarr_check(.false., 'real data into non-float array '//trim(arr%name))
        end if
    end subroutine zarr_write_chunk_3d_real

    ! ------------------------------------------------------------------ chunk writers (int input)

    subroutine zarr_write_chunk_1d_int(store, arr, cidx, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        integer,            intent(in) :: data(:)
        integer(int32), allocatable :: bi(:)
        integer :: c0, n0, i
        call zarr_check(trim(arr%dtype) == '<i4', 'int data into non-i4 array '//trim(arr%name))
        c0 = arr%chunks(1); n0 = size(data)
        allocate(bi(c0)); bi = int(nint(arr%fill), int32)
        do i = 1, n0; bi(i) = int(data(i), int32); end do
        call emit_i4(chunk_path(store, arr, cidx), bi, trim(arr%codec))
    end subroutine zarr_write_chunk_1d_int

    subroutine zarr_write_chunk_2d_int(store, arr, cidx, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: cidx(:)
        integer,            intent(in) :: data(:,:)
        integer(int32), allocatable :: bi(:)
        integer :: c0, c1, n0, n1, i0, i1
        call zarr_check(trim(arr%dtype) == '<i4', 'int data into non-i4 array '//trim(arr%name))
        c0 = arr%chunks(1); c1 = arr%chunks(2)
        n0 = size(data,1);  n1 = size(data,2)
        allocate(bi(c0*c1)); bi = int(nint(arr%fill), int32)
        do i0 = 1, n0
            do i1 = 1, n1
                bi((i0-1)*c1 + i1) = int(data(i0,i1), int32)
            end do
        end do
        call emit_i4(chunk_path(store, arr, cidx), bi, trim(arr%codec))
    end subroutine zarr_write_chunk_2d_int

    ! ------------------------------------------------------------------ chunk readers (RMW append)

    ! Read a chunk file into `raw` (nbytes decoded bytes), applying the codec inverse. Returns
    ! exist=.false. (and leaves raw unallocated) when the chunk file is absent — the caller then
    ! treats the chunk as all-fill (a fresh time-chunk). none => raw bytes verbatim; lz4 => strip the
    ! 4-byte numcodecs length header and LZ4_decompress_safe the block.
    subroutine read_chunk_bytes(path, nbytes, raw, exist, codec)
        character(len=*), intent(in)                            :: path, codec
        integer,          intent(in)                            :: nbytes
        integer(int8),    intent(out), allocatable, target      :: raw(:)
        logical,          intent(out)                           :: exist
        integer(int8),    allocatable, target :: filebytes(:)
        integer :: u, ios, fsz
        exist = .false.
        inquire(file=trim(path), exist=exist)
        if (.not. exist) return
        open(newunit=u, file=trim(path), status='old', action='read', &
             form='unformatted', access='stream', iostat=ios)
        call zarr_check(ios == 0, 'open(read chunk) '//trim(path))
        inquire(unit=u, size=fsz)
        allocate(filebytes(max(1,fsz)))
        if (fsz > 0) read(u) filebytes(1:fsz)
        close(u)
        select case (trim(codec))
        case ('none')
            call zarr_check(fsz == nbytes, 'read_chunk(none) size mismatch '//trim(path))
            call move_alloc(filebytes, raw)
        case ('lz4')
#ifdef HAVE_LZ4
            block
                integer(int32) :: declen
                integer(c_int) :: r
                call zarr_check(fsz >= 4, 'read_chunk(lz4) short header '//trim(path))
                declen = transfer(filebytes(1:4), declen)
                call zarr_check(int(declen) == nbytes, 'read_chunk(lz4) declen mismatch '//trim(path))
                allocate(raw(nbytes))
                r = c_lz4_decompress(c_loc(filebytes(5)), c_loc(raw(1)), &
                                     int(fsz - 4, c_int), int(nbytes, c_int))
                call zarr_check(int(r) == nbytes, 'LZ4_decompress_safe failed '//trim(path))
            end block
#else
            call zarr_check(.false., 'lz4 chunk read but HAVE_LZ4 off (relink with liblz4)')
#endif
        case default
            call zarr_check(.false., 'read_chunk unknown codec '//trim(codec))
        end select
    end subroutine read_chunk_bytes

    ! 1-D real chunk -> out(1:c0) (WP). Missing file => out = fill. Used for the time coord.
    ! The readers view raw as typed values (c_f_pointer), no transfer() copy (see emit_r8).
    subroutine zarr_read_chunk_1d_real(store, arr, cidx, out)
        type(t_zarr_store), intent(in)  :: store
        type(t_zarr_array), intent(in)  :: arr
        integer,            intent(in)  :: cidx(:)
        real(WP),           intent(out) :: out(:)
        integer(int8), allocatable, target :: raw(:)
        real(real64), pointer, contiguous :: b8(:)
        real(real32), pointer, contiguous :: b4(:)
        logical :: exist
        integer :: c0, nb
        c0 = arr%chunks(1)
        if (trim(arr%dtype) == '<f8') then; nb = 8*c0; else; nb = 4*c0; end if
        call read_chunk_bytes(chunk_path(store, arr, cidx), nb, raw, exist, trim(arr%codec))
        if (.not. exist) then; out(1:c0) = real(arr%fill, WP); return; end if
        if (trim(arr%dtype) == '<f8') then
            call c_f_pointer(c_loc(raw(1)), b8, [c0]); out(1:c0) = real(b8, WP)
        else
            call c_f_pointer(c_loc(raw(1)), b4, [c0]); out(1:c0) = real(b4, WP)
        end if
    end subroutine zarr_read_chunk_1d_real

    ! 2-D real chunk -> out(1:c0,1:c1) (WP), un-transposing C row-major. Missing => out = fill.
    subroutine zarr_read_chunk_2d_real(store, arr, cidx, out)
        type(t_zarr_store), intent(in)  :: store
        type(t_zarr_array), intent(in)  :: arr
        integer,            intent(in)  :: cidx(:)
        real(WP),           intent(out) :: out(:,:)
        integer(int8), allocatable, target :: raw(:)
        real(real64), pointer, contiguous :: b8(:)
        real(real32), pointer, contiguous :: b4(:)
        logical :: exist
        integer :: c0, c1, nb, i0, i1
        c0 = arr%chunks(1); c1 = arr%chunks(2)
        if (trim(arr%dtype) == '<f8') then; nb = 8*c0*c1; else; nb = 4*c0*c1; end if
        call read_chunk_bytes(chunk_path(store, arr, cidx), nb, raw, exist, trim(arr%codec))
        if (.not. exist) then; out(1:c0,1:c1) = real(arr%fill, WP); return; end if
        if (trim(arr%dtype) == '<f8') then
            call c_f_pointer(c_loc(raw(1)), b8, [c0*c1])
            do i0 = 1, c0; do i1 = 1, c1; out(i0,i1) = real(b8((i0-1)*c1 + i1), WP); end do; end do
        else
            call c_f_pointer(c_loc(raw(1)), b4, [c0*c1])
            do i0 = 1, c0; do i1 = 1, c1; out(i0,i1) = real(b4((i0-1)*c1 + i1), WP); end do; end do
        end if
    end subroutine zarr_read_chunk_2d_real

    ! 3-D real chunk -> out(1:c0,1:c1,1:c2) (WP), un-transposing C row-major. Missing => out = fill.
    subroutine zarr_read_chunk_3d_real(store, arr, cidx, out)
        type(t_zarr_store), intent(in)  :: store
        type(t_zarr_array), intent(in)  :: arr
        integer,            intent(in)  :: cidx(:)
        real(WP),           intent(out) :: out(:,:,:)
        integer(int8), allocatable, target :: raw(:)
        real(real64), pointer, contiguous :: b8(:)
        real(real32), pointer, contiguous :: b4(:)
        logical :: exist
        integer :: c0, c1, c2, nb, i0, i1, i2
        c0 = arr%chunks(1); c1 = arr%chunks(2); c2 = arr%chunks(3)
        if (trim(arr%dtype) == '<f8') then; nb = 8*c0*c1*c2; else; nb = 4*c0*c1*c2; end if
        call read_chunk_bytes(chunk_path(store, arr, cidx), nb, raw, exist, trim(arr%codec))
        if (.not. exist) then; out(1:c0,1:c1,1:c2) = real(arr%fill, WP); return; end if
        if (trim(arr%dtype) == '<f8') then
            call c_f_pointer(c_loc(raw(1)), b8, [c0*c1*c2])
            do i0 = 1, c0; do i1 = 1, c1; do i2 = 1, c2
                out(i0,i1,i2) = real(b8((i0-1)*c1*c2 + (i1-1)*c2 + i2), WP)
            end do; end do; end do
        else
            call c_f_pointer(c_loc(raw(1)), b4, [c0*c1*c2])
            do i0 = 1, c0; do i1 = 1, c1; do i2 = 1, c2
                out(i0,i1,i2) = real(b4((i0-1)*c1*c2 + (i1-1)*c2 + i2), WP)
            end do; end do; end do
        end if
    end subroutine zarr_read_chunk_3d_real

    ! ------------------------------------------------------------------ whole-array writers
    ! Convenience for the 1-rank / single-writer path: write every chunk of a full in-memory array.
    ! Multi-rank writers call zarr_write_chunk directly (each writes only its own chunks).

    subroutine zarr_write_whole_1d_real(store, arr, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(real64),       intent(in) :: data(:)
        integer :: s0, c0, ic0, a0, n0
        s0 = arr%dims(1); c0 = arr%chunks(1)
        call zarr_check(size(data) == s0, 'write_whole_1d size mismatch '//trim(arr%name))
        do ic0 = 0, (s0 + c0 - 1)/c0 - 1
            a0 = ic0*c0 + 1; n0 = min(c0, s0 - ic0*c0)
            call zarr_write_chunk_1d_real(store, arr, [ic0], data(a0:a0+n0-1))
        end do
    end subroutine zarr_write_whole_1d_real

    subroutine zarr_write_whole_2d_real(store, arr, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        real(real64),       intent(in) :: data(:,:)
        integer :: s0, s1, c0, c1, ic0, ic1, a0, b0, n0, n1
        s0 = arr%dims(1); s1 = arr%dims(2); c0 = arr%chunks(1); c1 = arr%chunks(2)
        call zarr_check(size(data,1) == s0 .and. size(data,2) == s1, &
                        'write_whole_2d size mismatch '//trim(arr%name))
        do ic0 = 0, (s0 + c0 - 1)/c0 - 1
            a0 = ic0*c0 + 1; n0 = min(c0, s0 - ic0*c0)
            do ic1 = 0, (s1 + c1 - 1)/c1 - 1
                b0 = ic1*c1 + 1; n1 = min(c1, s1 - ic1*c1)
                call zarr_write_chunk_2d_real(store, arr, [ic0,ic1], data(a0:a0+n0-1, b0:b0+n1-1))
            end do
        end do
    end subroutine zarr_write_whole_2d_real

    subroutine zarr_write_whole_1d_int(store, arr, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: data(:)
        integer :: s0, c0, ic0, a0, n0
        s0 = arr%dims(1); c0 = arr%chunks(1)
        call zarr_check(size(data) == s0, 'write_whole_1d size mismatch '//trim(arr%name))
        do ic0 = 0, (s0 + c0 - 1)/c0 - 1
            a0 = ic0*c0 + 1; n0 = min(c0, s0 - ic0*c0)
            call zarr_write_chunk_1d_int(store, arr, [ic0], data(a0:a0+n0-1))
        end do
    end subroutine zarr_write_whole_1d_int

    subroutine zarr_write_whole_2d_int(store, arr, data)
        type(t_zarr_store), intent(in) :: store
        type(t_zarr_array), intent(in) :: arr
        integer,            intent(in) :: data(:,:)
        integer :: s0, s1, c0, c1, ic0, ic1, a0, b0, n0, n1
        s0 = arr%dims(1); s1 = arr%dims(2); c0 = arr%chunks(1); c1 = arr%chunks(2)
        call zarr_check(size(data,1) == s0 .and. size(data,2) == s1, &
                        'write_whole_2d size mismatch '//trim(arr%name))
        do ic0 = 0, (s0 + c0 - 1)/c0 - 1
            a0 = ic0*c0 + 1; n0 = min(c0, s0 - ic0*c0)
            do ic1 = 0, (s1 + c1 - 1)/c1 - 1
                b0 = ic1*c1 + 1; n1 = min(c1, s1 - ic1*c1)
                call zarr_write_chunk_2d_int(store, arr, [ic0,ic1], data(a0:a0+n0-1, b0:b0+n1-1))
            end do
        end do
    end subroutine zarr_write_whole_2d_int

end module mod_io_zarr
