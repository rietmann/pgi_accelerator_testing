home=`pwd`
for dir in `ls --color=never -d */`
do
	echo $dir

	for thd in 1 2 4 6 12 24
	do
		thddir=${home}/${dir}/thd${thd}
		echo ${thddir}

		mkdir ${thddir}
		cp ${home}/${dir}/bench ${thddir}

		# create the autotune.sh script
		# requires that a parameter file, "autotune.params", is located in the directory
		file=${thddir}/autotune.sh
		echo java -jar ${home}/patus.jar qsubautotune ${thddir}/bench ${home}/submit-template.job ${thd} `cat ${home}/${dir}/autotune.params` \> ${thddir}/result.txt > ${file}

		# execute the script
		chmod +x ${thddir}/autotune.sh
		cd ${thddir}
		#nohup ./autotune.sh &
		./autotune.sh

		#exit
	done
done
