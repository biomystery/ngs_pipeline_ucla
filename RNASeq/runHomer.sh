rg=10000
dir="./motif_$rg"
mkdir $dir
for ((i=1; i<=10;i++))
do
  findMotifs.pl clust"$i".txt mouse $dir/clust"$i"  -start -$rg -end 200 -p 15 -S 15 > $dir/out"$i".txt &
done
cp runHomer.sh $dir/ 