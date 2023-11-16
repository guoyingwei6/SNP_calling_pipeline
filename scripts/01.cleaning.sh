${sm}=$1
fastp -i ${sm}_1.fq.gz -I ${sm}_2.fq.gz -o ${sm}_1.clean.fq.gz -O ${sm}_2.clean.fq.gz -h ${sm}.html -j ${sm}.json -w 8