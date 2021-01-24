#!/bin/sh

echo splitting 538 frames...
convert ../galton_board-538.gif -coalesce a-%04d.gif

echo splitting economist frames...
convert ../galton_board-the_economist.gif -coalesce b-%04d.gif

echo merging frames...
for f in a-*.gif; do
    echo $f
    convert $f ${f/a/b} +append $f
done

echo making gif...
convert -loop 0 -delay 3 a-*.gif galton_board-side_by_side.gif

echo optimizing gif...
mogrify -layers 'optimize' -fuzz 3% galton_board-side_by_side.gif

rm a*.gif b*.gif

