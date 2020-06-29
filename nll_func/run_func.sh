echo "============"
echo "Script to run NLLoc_func_test program"
echo "uses phase data from NonLinLoc Sample-Location Tutorial for Local and Regional Events"
echo
echo "IMPORTANT:  Requires:"
echo "   1. NonLinLoc - the command \"NLLoc_func_test\" must be on your path"
echo "   2. Java - the command \"java\" must be on your path"
echo "   3. SeismicityViewer must be installed and on your java classpath - see: http://alomax.net/seismicity"
echo
echo "See: NLLoc_func_test.c for more information on using NLLoc through a function call."
echo "============"
rm out/*
NLLoc_func_test run/nlloc_func_test.in obs/vinti_1.obs obs/vinti_2.obs obs/vinti_3.obs obs/vinti_4.obs obs/vinti_5.obs
java -Xmx768m net.alomax.seismicity.Seismicity out/*.hyp
