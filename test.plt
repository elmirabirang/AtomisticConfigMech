set terminal svg enhanced size 1000 1000 fname "Times" fsize 36
set output "plot.svg"
set title "A simple plot of x^2 vs. x"
set xlabel "x"
set ylabel "y"
set zlabel "z"

#splot "./test.data" using 3:1:2 title ""

splot "< awk '{if($4 == \"red\") print}' test.data" u 3:1:2 t "" w p pt 7, \
      "< awk '{if($4 == \"blue\") print}' test.data" u 3:1:2 t "" w p pt 7,\
      "< awk '{if($4 == \"yellow\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"green\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"pink\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"violet\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"black\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"orange\") print}' test.data" u 3:1:2 t "" w p pt 7 ,\
      "< awk '{if($4 == \"cyan\") print}' test.data" u 3:1:2 t "" w p pt 7


