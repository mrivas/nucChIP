for file in /var/www/public/nucChIP/files/fragmentDist/exons/*rst; do
	ln -f -s $file .
done

for file in /var/www/public/nucChIP/files/fragmentDist/tss/*rst; do
	ln -f -s $file .
done
