RE="125 140 150 175 200 250 300 400 500 600"
TAU="80 90 100 105"
NCRIT="5 6 7 8 9"
FOIL=dai1336a
for r in $RE
do
    for t in $TAU
    do
        for n in $NCRIT
        do
            ./genpolar.sh $FOIL.$t $r $n
        done
    done
done
