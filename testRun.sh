spray-and-pray.py -g toy_dataset/test_dataset.fasta -ref toy_dataset/test_database.faa --makedb --bin -out testrun --meta -t 4 --test
md5sum-check.py -i toy_dataset -o testrun
