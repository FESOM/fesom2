# XIOS integration sketch for FESOM-2 output

Goal: route CMIP7 output (especially heavy daily 3D r8 fields like `utemp`, `vtemp`) through
XIOS async server ranks, replacing the single-root serial `nf90_put_var` path in
`src/io_meandata.F90` that caps daily cost at ~6 s/day.

Assumptions:
- XIOS 2.x source is already in-tree at `../xios`. Build + link to produce `libxios.a`.
- OpenMPI 4.1 on Levante provides `MPI_THREAD_MULTIPLE` (kept for restart/other code paths; XIOS itself doesn't need it).
- Run launch uses a hybrid layout: N FESOM compute ranks + M XIOS server ranks in the same
  MPI job (OASIS/launcher already does this kind of split for IFS/XIOS; same mechanism reused).

Not in scope here: XIOS *client-server* parallel output, just a standalone FESOM→XIOS path.
Coupled behaviour (AWI-ESM3 with OIFS) stays on existing MULTIO path until we know XIOS
works standalone.

---

## 1. Build

### 1.0. Rebuild via esm_tools

From the parent dir:
```
cd /work/ab0246/a270092/model_codes
esm_master comp-awiesm3-develop/fesom      # incremental recompile (use this day-to-day)
# first-time only: esm_master install-awiesm3-develop/fesom
```
The `with_xios` toggle is awiesm3-local (`configs/setups/awiesm3/awiesm3.yaml`, default `true`), which injects `-DFESOM_WITH_XIOS=ON` + `export XIOS_ROOT=${model_dir}/../xios` into the fesom `comp_command` via a `choose_with_xios` block. Standalone `components/fesom/fesom-2.7.yaml` is unchanged.

### 1a. XIOS
```
cd xios
./make_xios --arch GCC_LEVANTE --job 16         # or whatever arch file matches Levante
# produces lib/libxios.a and bin/xios_server.exe
```

### 1b. FESOM CMake changes (`src/CMakeLists.txt`)
Add option + link:
```cmake
option(FESOM_WITH_XIOS "Enable XIOS output" OFF)
if(FESOM_WITH_XIOS)
    find_path(XIOS_INCLUDE_DIR NAMES xios.mod HINTS $ENV{XIOS_ROOT}/inc)
    find_library(XIOS_LIB NAMES xios HINTS $ENV{XIOS_ROOT}/lib)
    target_include_directories(fesom PRIVATE ${XIOS_INCLUDE_DIR})
    target_link_libraries(fesom PRIVATE ${XIOS_LIB} stdc++)
    target_compile_definitions(fesom PRIVATE __XIOS)
endif()
```

---

## 2. New module: `src/io_xios.F90`

Skeleton (names approximate; adapt to FESOM style):

```fortran
module io_xios_module
#if defined(__XIOS)
  use xios
  use mod_mesh
  use mod_partit
  use g_clock, only: timenew
  implicit none
  private

  public :: xios_fesom_init, xios_fesom_close
  public :: xios_fesom_update_calendar, xios_fesom_send_field_2d, xios_fesom_send_field_3d

  type(xios_context)    :: ctx
  integer, save         :: xios_comm      ! FESOM client comm returned by XIOS
  logical, save         :: xios_on = .false.

contains

  subroutine xios_fesom_init(mesh, partit, parent_comm, out_comm)
    type(t_mesh),   intent(in)    :: mesh
    type(t_partit), intent(in)    :: partit
    integer,        intent(in)    :: parent_comm
    integer,        intent(out)   :: out_comm
    integer                       :: i, ierr
    real(8), allocatable          :: lon(:), lat(:)
    integer, allocatable          :: i_index(:)
    real(8), allocatable          :: axis_values(:)

    call xios_initialize("fesom", return_comm=out_comm)   ! splits off XIOS servers
    xios_comm = out_comm

    call xios_context_initialize("fesom_ctx", xios_comm)
    call xios_get_handle("fesom_ctx", ctx)
    call xios_set_current_context(ctx)

    ! --- unstructured node-based domain ---
    allocate(lon(partit%myDim_nod2D), lat(partit%myDim_nod2D))
    allocate(i_index(partit%myDim_nod2D))
    do i = 1, partit%myDim_nod2D
       lon(i)     = mesh%coord_nod2D(1, i) * (180d0 / pi)
       lat(i)     = mesh%coord_nod2D(2, i) * (180d0 / pi)
       i_index(i) = partit%myList_nod2D(i) - 1   ! 0-based
    end do

    call xios_set_domain_attr("nodes",                 &
         type = "unstructured",                         &
         ni_glo  = mesh%nod2D,                          &
         ni      = partit%myDim_nod2D,                  &
         ibegin  = 0,                                   &
         i_index = i_index,                             &
         lonvalue_1d = lon, latvalue_1d = lat,          &
         nvertex = 1)                                   ! use 3 + bounds for CF UGRID

    ! --- element-based domain (for u, v, vorticity, ...) ---
    ! ...analogous with mesh%elem2D, myList_elem2D, coord_nod2D(:, elem2D_nodes(:,:)) for bounds

    ! --- vertical axis ---
    allocate(axis_values(mesh%nl-1))
    do i = 1, mesh%nl-1
       axis_values(i) = 0.5d0 * (mesh%zbar(i) + mesh%zbar(i+1))
    end do
    call xios_set_axis_attr("nz", n_glo = mesh%nl-1, value = axis_values)

    ! --- calendar (must match FESOM's clock) ---
    call xios_define_calendar(type="D360" /* or Gregorian */, &
         start_date = xios_date(yearstart, 1, 1, 0, 0, 0),    &
         time_origin = xios_date(yearstart, 1, 1, 0, 0, 0))

    call xios_close_context_definition()
    xios_on = .true.

    deallocate(lon, lat, i_index, axis_values)
  end subroutine xios_fesom_init


  subroutine xios_fesom_update_calendar(step)
    integer, intent(in) :: step
    if (.not. xios_on) return
    call xios_set_current_context(ctx)
    call xios_update_calendar(step)
  end subroutine


  subroutine xios_fesom_send_field_2d(name, local_vals)
    character(len=*), intent(in) :: name
    real(8),          intent(in) :: local_vals(:)     ! shape (myDim_nod2D)
    if (.not. xios_on) return
    call xios_send_field(trim(name), local_vals)
  end subroutine


  subroutine xios_fesom_send_field_3d(name, local_vals)
    character(len=*), intent(in) :: name
    real(8),          intent(in) :: local_vals(:,:)   ! shape (nz, myDim_nod2D) or swapped
    if (.not. xios_on) return
    call xios_send_field(trim(name), local_vals)
  end subroutine


  subroutine xios_fesom_close()
    if (.not. xios_on) return
    call xios_context_finalize()
    call xios_finalize()
  end subroutine

#endif
end module io_xios_module
```

Notes:
- Shape of `local_vals` must match the grid/axis order declared in XML. XIOS is fussy
  here — wrong ordering gives a silent transpose at the output.
- For element-based fields (if any), define a second domain `"elements"` and a second
  grid `"grid_3d_elem"`. Many CMIP7 ocean variables are node-based so this is secondary.

---

## 3. Hooks in existing FESOM code

### 3a. Startup (`src/fesom_module.F90`)

After MPI init, OASIS init, mesh read, and partitioning — but before the main time loop:

```fortran
#if defined(__XIOS)
   call xios_fesom_init(f%mesh, f%partit, f%partit%MPI_COMM_FESOM, xios_client_comm)
   ! replace f%partit%MPI_COMM_FESOM with xios_client_comm from here on
   ! (all FESOM MPI collectives must use the shrunk communicator)
#endif
```

**Critical subtlety**: `xios_initialize` splits the parent communicator into
(FESOM clients) + (XIOS servers). From that point on, FESOM's own communicator must be
the returned client comm, not MPI_COMM_WORLD. This is the single most error-prone part
of the integration — miss it and half the ranks hang inside `MPI_Allreduce`.

For AWI-ESM3 coupled runs, OASIS already does a WORLD→model split. XIOS needs to split
FESOM's model-local communicator further. Check that OASIS3-MCT's `oasis_get_localcomm`
is called *before* `xios_initialize`, and that the latter splits that local comm.

### 3b. Timestep loop (`src/fesom_module.F90`, inside `step(...)`)

```fortran
#if defined(__XIOS)
   call xios_fesom_update_calendar(istep)
#endif
```

Place once per step, at the top of `step()` or right before the output section.

### 3c. Output path (`src/io_meandata.F90`, inside `output()`)

Option A — **parallel pilot**: keep existing path, *additionally* send to XIOS for a
chosen subset of streams. Controlled by a namelist flag per stream or globally.

Inside the per-stream loop, after the accumulated `local_values_r8_copy` is computed
(around line 2500), add:

```fortran
#if defined(__XIOS)
   if (entry%use_xios) then
      if (entry%accuracy == i_real8) then
         if (size(entry%local_values_r8_copy,1) == 1) then
            call xios_fesom_send_field_2d(entry%name, entry%local_values_r8_copy(1,:))
         else
            call xios_fesom_send_field_3d(entry%name, entry%local_values_r8_copy)
         end if
      else
         ! analogous for r4
      end if
      cycle     ! skip the serial nf90_put_var path
   end if
#endif
```

Option B — **cutover**: once parallel pilot proves correctness + speed, remove the
existing nf90 write path for those streams entirely.

### 3d. Shutdown (`src/fesom_module.F90`)

```fortran
#if defined(__XIOS)
   call xios_fesom_close()
#endif
```

Before `MPI_Finalize`.

---

## 4. XML: `iodef.xml`

Minimal example controlling frequency, averaging, and filename:

```xml
<?xml version="1.0"?>
<simulation>
  <context id="fesom_ctx" calendar_type="Gregorian" start_date="1900-01-01 00:00:00">

    <axis_definition>
      <axis id="nz" name="depth" unit="m"/>
    </axis_definition>

    <domain_definition>
      <domain id="nodes" />
      <domain id="elements" />
    </domain_definition>

    <grid_definition>
      <grid id="grid_2d_nod"> <domain domain_ref="nodes"/>       </grid>
      <grid id="grid_3d_nod"> <domain domain_ref="nodes"/> <axis axis_ref="nz"/> </grid>
    </grid_definition>

    <field_definition prec="8" operation="average" freq_op="30mi">
      <field id="temp"   name="thetao" grid_ref="grid_3d_nod" unit="degC"/>
      <field id="salt"   name="so"     grid_ref="grid_3d_nod" unit="psu"/>
      <field id="ssh"    name="zos"    grid_ref="grid_2d_nod" unit="m"/>
      <field id="utemp"  name="utemp"  grid_ref="grid_3d_nod" unit="degC m s-1"/>
      <field id="vtemp"  name="vtemp"  grid_ref="grid_3d_nod" unit="degC m s-1"/>
      <!-- … the rest … -->
    </field_definition>

    <file_definition type="one_file" output_freq="1d" sync_freq="10d">
      <file id="daily_3d" name="fesom_daily_3d">
        <field field_ref="temp"/>
        <field field_ref="salt"/>
        <field field_ref="utemp"/>
        <field field_ref="vtemp"/>
      </file>
      <file id="daily_2d" name="fesom_daily_2d">
        <field field_ref="ssh"/>
        <!-- … -->
      </file>
    </file_definition>

  </context>
</simulation>
```

`type="one_file"` → XIOS servers write a single NetCDF per output file (most CMIP-like).
`type="multiple_file"` splits per-server for peak throughput; we'd combine in
postproc. Start with `one_file` and benchmark.

---

## 5. Launch layout

Batch script `srun` layout with server ranks, e.g. on 27 compute nodes / 768 FESOM ranks:

```
# 24 FESOM nodes + 3 XIOS server nodes (~3% of resource for I/O)
srun --multi-prog mpmd.conf
```

`mpmd.conf`:
```
0-767       ./fesom.x
768-959     ./xios_server.exe
```

Start with a 32:1 ratio (24 servers for 768 compute) and tune.

---

## 6. Sanity milestones

1. XIOS build + link clean, `FESOM_WITH_XIOS=ON` compiles, run with `xios_on=false` still behaves identically.
2. `xios_fesom_init` + domain definition completes without hanging (no send, just close after init).
3. Send **one** 2D node field (e.g. `ssh`) and get a valid UGRID NetCDF output.
4. Compare `cdo diff ssh_xios.nc ssh_old.nc` — bitwise or at least <1e-12 agreement.
5. Add 3D node fields; benchmark `utemp` daily write cost.
6. Roll out remaining streams, then delete the unused streams from the old `io_meandata` path.

---

## 7. Likely pitfalls (known from NEMO / IPSL-CM experience)

- **Communicator mistake** (section 3a). Symptom: some ranks hang in the first collective after init.
- **Data shape mismatch**. XIOS transposes silently if the shape you send doesn't match the grid definition. Symptom: spatial nonsense in the output, no error.
- **Calendar drift** between FESOM's `g_clock` and XIOS's. Symptom: output timestamps off by 1 step.
- **Non-contiguous `myList_nod2D`**. XIOS expects a 0-based `i_index` array; FESOM's list is 1-based, so subtract 1 on copy.
- **Halo nodes**. Only send `1:myDim_nod2D` (owned), not the full array including halos. The ghost values are stale.
- **XIOS version/arch file**. `arch/arch-GCC_LEVANTE.*` may or may not exist; may need to create from `arch-GCC_MUNICH` or similar.

---

## 8. Path forward

Suggested phasing:
1. Week 1: build + wiring (sections 1-3), get `xios_on=false` path stable, get UGRID output of `ssh` matching old output.
2. Week 2: add `temp`, `salt`, `utemp`, `vtemp`. Benchmark daily cost. Decide ratio of XIOS servers.
3. Week 3: port remaining streams, write CMIP7 `iodef.xml` subset, regression-test year-long run.
4. Week 4: integrate with AWI-ESM3 coupled config (OASIS + MULTIO + XIOS coexistence).

Fallback at any stage: `FESOM_WITH_XIOS=OFF` or per-stream `use_xios=.false.` → existing path.
