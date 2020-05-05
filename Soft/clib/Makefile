bs2px.so: bs2px.c
	gcc -shared -fpic -fopenmp -o bs2px.so bs2px.c
map_making.so: map_making.c
	gcc -shared -fpic -o map_making.so map_making.c
bm_convolve.so: bm_convolve.c
	gcc -shared -fpic -o bm_convolve.so bm_convolve.c
grd_sub.so: grd_sub.c
	gcc -shared -fpic -fopenmp -o grd_sub.so grd_sub.c
map_accumulation.so: map_accumulation.c
	gcc -shared -fpic -fopenmp -o map_accumulation.so map_accumulation.c


.PHONY:clean
clean:$(objects)
	rm $(objects)
