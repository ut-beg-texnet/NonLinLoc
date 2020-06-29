



echo "---------------------------------"
geo_test -49.716797 -114.931641  -38.0567 -61.9795
echo "should be:"
echo "distance: 4334.965466151394 (km)"
echo "azimuth: 93.33977502156569"
echo

echo "---------------------------------"
geo_test 39.5654 -167.241  33.0033 136.779
echo "should be:"
echo "distance: 4986.077239913946 (km)"
echo "azimuth: 279.8776594786318"
echo

echo "---------------------------------"
geo_test 0 0 0 0
echo "should be:"
echo "distance: 0.0 (km)"
echo "azimuth: 0.0"
echo

echo "---------------------------------"
geo_test 0 0 0 90
echo "should be:"
echo "distance: 10000.0 (km)"
echo "azimuth: 90.0"
echo

