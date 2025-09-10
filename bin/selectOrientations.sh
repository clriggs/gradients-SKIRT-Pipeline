#!/bin/bash
cptmarvel_halos=(1 2 3 5 6 7 10 11 13)
elektra_halos=(1 2 3 4 5 9 10 11 12 17 36)
storm_halos=(1 2 3 4 5 6 7 8 10 11 12 14 15 23 31 44)
rogue_halos=(1 3 7 8 10 11 12 15 16 17 28 31 37)
h148_halos=(2 3 4 6 7 11 12 13 15 20 23 27 28 33 34 38 65 86 114)
h229_halos=(2 3 6 14 18 22 49 92)
h242_halos=(8 10 21 30 34 38)
h329_halos=(7 29 115 117)

mdcjl_path="/data/riggs/gradients/marvel_dcjl_sim_dict.pickle"

declare -A mdcjl_dict

mdcjl_dict[cptmarvel]=cptmarvel_halos[@]
mdcjl_dict[elektra]=elektra_halos[@]
mdcjl_dict[storm]=storm_halos[@]
mdcjl_dict[rogue]=rogue_halos[@]
mdcjl_dict[h148]=h148_halos[@]
mdcjl_dict[h229]=h229_halos[@]
mdcjl_dict[h242]=h242_halos[@]
mdcjl_dict[h329]=h329_halos[@]

for i in "${!mdcjl_dict[@]}"; do
    echo "on sim: $i"
    for j in "${!mdcjl_dict[$i]}"; do
        echo "on halo: $j"
        python selectOrientations.py --sim=$i --sim_dict_path=$mdcjl_path --halo=$j
    done
done