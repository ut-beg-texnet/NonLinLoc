

rm -r nlloc_sample_test
mkdir nlloc_sample_test
cp -pr nlloc_sample/* nlloc_sample_test
cd nlloc_sample_test
pwd
setGMT4
./run_nlloc_sample.sh
echo "Following should indicate no difference:"
COMMAND="diff original_output/alaska.hyp ../nlloc_sample_test_frozen_20181208/original_output/alaska.hyp"
echo ${COMMAND}
${COMMAND}
cd ..

