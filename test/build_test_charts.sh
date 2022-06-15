#!/bin/bash

for distance in 1500 2000 2500; do
  for size in a0 a1 a2 a3 a4; do
    ./mtf_generate_test_chart -t perspective -s $size -d $distance -o temp.svg
    inkscape --export-pdf=perspective_${size}_distance_${distance}.pdf temp.svg
  done
done
rm temp.svg

for size in a0 a1 a2 a3 a4; do
  ./mtf_generate_test_chart -t grid -s $size -o temp.svg
  inkscape --export-pdf=grid_${size}.pdf temp.svg
done
rm temp.svg

for size in a0 a1 a2 a3 a4; do
  ./mtf_generate_test_chart -t lensgrid -s $size -o temp.svg
  inkscape --export-pdf=lensgrid_${size}.pdf temp.svg
done
rm temp.svg

size=a3
./mtf_generate_test_chart -t mfperspective -o temp.svg
inkscape --export-pdf=mfperspective_${size}.pdf temp.svg
rm temp.svg

./mtf_generate_test_chart -t focus  -o temp.svg
inkscape --export-pdf=focus_${size}.pdf temp.svg
rm temp.svg
