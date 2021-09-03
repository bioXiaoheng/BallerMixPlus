# script to pick out the alignment terms relevant to a certain chromosome

input=$1
chrom=$2

zcat $input | awk -v c=$chrom 'BEGIN{ counter = 4 }{
	if (!($0 ~ /^#/)){
		counter = counter - 1
		if ( counter == 3){
			if ($2 == "chr"c && $5 ~ /^(chr)?([0-9]+|[X|x]|[Y|y])[a-z|A-Z]*$/  ){
				ok = 1
				print $0
			}else{
				ok = 0
			}
		}else if (counter == 0){
			counter = 4
			if(ok == 1){
				print $0
			}
		}else{
			if(ok == 1){
				print $0
			}
		}
	}
}'