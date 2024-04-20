Wall_Collisions=0

D_Manual=1
D_Perp=0 
D_Par=0.00001
D_Rot=0 

Skip=1000

mkdir -p runs

cd runs
for seed in {1..800}
do 

mkdir -p run${seed}

cd run${seed}

mkdir -p plots
mkdir -p data

cp ../../template.inp input.inp
cp ../../filaments . 
cp ../../reader.py .

sed -i -e "s/VAL_SEED/${seed}/g" input.inp

sed -i -e "s/VAL_COLLISIONS/${Wall_Collisions}/g" input.inp

sed -i -e "s/VAL_DMANUAL/${D_Manual}/g" input.inp
sed -i -e "s/VAL_DPERP/${D_Perp}/g" input.inp
sed -i -e "s/VAL_DPAR/${D_Par}/g" input.inp
sed -i -e "s/VAL_DROT/${D_Rot}/g" input.inp

sed -i -e "s/VAL_SKIP/${Skip}/g" input.inp





./filaments input.inp
# python3 reader.py

echo '-----------------------------------'
echo

cd ../
done
cd ../  