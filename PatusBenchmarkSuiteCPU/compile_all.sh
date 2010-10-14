home=`pwd`
for f in `ls --color=never -d */`
do
	echo '-->' $f '...'
	cd ${home}/${f}
	make $1
done
