#!/bin/bash

new_folder="static_pattern/"
new_folder="uniform/"

theta_arr="theta1/ theta1.5/ theta2.5/"
seed_arr="seed_half/ seed_single/"

mkdir -p "$new_folder"
echo Folder "$new_folder" is created.
for seed in $seed_arr
do
 mkdir "${new_folder}${seed}"
 echo Folder "${new_folder}${seed}" is created
 for theta in $theta_arr
 do
  mkdir "$new_folder${seed}${theta}"
  echo Folder "$new_folder${seed}${theta}" is created
 done
done
