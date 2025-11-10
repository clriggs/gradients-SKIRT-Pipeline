#!/bin/bash
r431_halos=(1)
r431_sidm_halos=(1)
r442_halos=(1)
r468_halos=(1)
r488_halos=(1)
r489_halos=(1)
r492_halos=(1)
r492_sidm_halos=(1)
r502_halos=(1)
r515_halos=(1)
r523_halos=(1)
r544_halos=(1)
r552_halos=(1)
r555_halos=(1)
r556_halos=(1)
r563_halos=(1)
r568_halos=(1)
r569_halos=(1)
r571_halos=(1)
r597_halos=(1)
r613_halos=(1)
r614_halos=(1)
r615_halos=(1)
r618_halos=(1)
r634_halos=(1)
r642_halos=(1)
r656_halos=(1)
r707_halos=(1)
r716_halos=(1)
r718_halos=(1)
r719_halos=(1)
r741_halos=(1)
r749_halos=(1)
r753_halos=(1)
r761_halos=(1)
r850_halos=(1)
r852_halos=(1)
r886_halos=(1)
r916_halos=(1)
r918_halos=(1)
r968_halos=(3)
r977_halos=(1)
r1023_halos=(3)

rom_path="/data/riggs/gradients-SKIRT-Pipeline/resources/romulus25_sim_dict.pickle"

declare -A rom_dict

rom_dict[r431]=r431_halos[@]
rom_dict[r431_sidm]=r431_sidm_halos[@]
rom_dict[r442]=r442_halos[@]
rom_dict[r468]=r468_halos[@]
rom_dict[r488]=r488_halos[@]
rom_dict[r489]=r489_halos[@]
rom_dict[r492]=r492_halos[@]
rom_dict[r492_sidm]=r492_sidm_halos[@]
rom_dict[r502]=r502_halos[@]
rom_dict[r515]=r515_halos[@]
rom_dict[r523]=r523_halos[@]
rom_dict[r544]=r544_halos[@]
rom_dict[r552]=r552_halos[@]
rom_dict[r555]=r555_halos[@]
rom_dict[r556]=r556_halos[@]
rom_dict[r563]=r563_halos[@]
rom_dict[r568]=r568_halos[@]
rom_dict[r569]=r569_halos[@]
rom_dict[r571]=r571_halos[@]
rom_dict[r597]=r597_halos[@]
rom_dict[r613]=r613_halos[@]
rom_dict[r614]=r614_halos[@]
rom_dict[r615]=r615_halos[@]
rom_dict[r618]=r618_halos[@]
rom_dict[r634]=r634_halos[@]
rom_dict[r642]=r642_halos[@]
rom_dict[r656]=r656_halos[@]
rom_dict[r707]=r707_halos[@]
rom_dict[r716]=r716_halos[@]
rom_dict[r718]=r718_halos[@]
rom_dict[r719]=r719_halos[@]
rom_dict[r741]=r741_halos[@]
rom_dict[r749]=r749_halos[@]
rom_dict[r753]=r753_halos[@]
rom_dict[r761]=r761_halos[@]
rom_dict[r850]=r850_halos[@]
rom_dict[r852]=r852_halos[@]
rom_dict[r886]=r886_halos[@]
rom_dict[r916]=r916_halos[@]
rom_dict[r918]=r918_halos[@]
rom_dict[r968]=r968_halos[@]
rom_dict[r977]=r977_halos[@]
rom_dict[r1023]=r1023_halos[@]

for i in "${!rom_dict[@]}"; do
    echo "on sim: $i"
    for j in "${!rom_dict[$i]}"; do
        echo "on halo: $j"
        python sampleOrientations.py --num=10 --sim=$i --sim_dict_path=$rom_path --halo=$j --SF=True --tauClear=3 --FoV=74240 --distance=50
    done
done