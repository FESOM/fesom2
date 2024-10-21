# Download icepack columphysics code
# from Lorenzo's fork on github and
# switch to the appropriate branch

DIR="./Icepack"
if [ ! -d "$DIR" ]; then
   git clone https://github.com/lzampier/Icepack.git
   cd $DIR
   git checkout icepack_fesom2
fi
