home=`pwd`
for dir in `ls --color=never -d */`
do
	echo $dir
	rm ${home}/${dir}/*.txt
	rm ${home}/${dir}/*.err
	rm ${home}/${dir}/*.job

	for thd in 1 2 4 6 12 24
	do
		thddir=${home}/${dir}/thd${thd}
		echo ${thddir}

		rm -rf ${thddir}
	done
done
