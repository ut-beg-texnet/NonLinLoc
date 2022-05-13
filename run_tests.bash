

rm -r nlloc_sample_test
mkdir nlloc_sample_test
cp -pr nlloc_sample/* nlloc_sample_test
cd nlloc_sample_test
pwd
setGMT4
./run_nlloc_sample.sh

cp -pr loc/* original_output/

echo ""
echo "----------------------------"
echo "Verify run dates:"
COMMAND="grep SIGNATURE original_output/alaska.hyp ../nlloc_sample_test_frozen_20181208/original_output/alaska.hyp"
echo ${COMMAND}
${COMMAND}

echo ""
echo "----------------------------"
echo "Following should indicate no significant differences in results:"
COMMAND="diff ../nlloc_sample_test_frozen_20220513/original_output/alaska.hyp original_output/alaska.hyp"
echo ""
echo "${COMMAND} | grep GEOGRAPHIC"
${COMMAND} | grep GEOGRAPHIC
echo ""
echo "${COMMAND} | grep QML_OriginUncertainty"
${COMMAND} | grep QML_OriginUncertainty

cd ..

