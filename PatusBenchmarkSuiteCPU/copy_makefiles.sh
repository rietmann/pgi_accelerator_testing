home=`pwd`
for f in `ls --color=never -d */`
do
	echo Copying ${home}/Makefile to ${home}/${f}/Makefile...
	cp ${home}/Makefile ${home}/${f}/Makefile 
done 
