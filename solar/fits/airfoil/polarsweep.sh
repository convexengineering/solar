RE="500 1000 1500 2000 2500 3000"
#RE="1500"
#TAU="80 90 100 110 120 130 140 150 160 170 180 190"
#TAU="200 210 220 230"
#TAU="120"
FOIL=dai1336
for r in $RE
do
    for t in $TAU
    do
        ./genpolar.sh $FOIL $t $r
    done
done
