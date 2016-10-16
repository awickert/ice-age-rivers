# Example
for name in ICE6G ICE5G ICE3G ANU
do
  convert -delay 50 -reverse map_$name*.png $name.mov
done

name=G12
convert -delay 10 -reverse map_$name*.png $name.mov

