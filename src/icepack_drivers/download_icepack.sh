# Download icepack columphysics code
# from Frank's fork on github and
# switch to the appropriate branch

DIR="./Icepack"
if [ ! -d "$DIR" ]; then
   git clone https://github.com/fkauker/Icepack.git
   cd $DIR
   git checkout icepack1.5.0_fesom2
fi
