
cat << END > alaska.xy
> GMT_LONLAT
END
pscoast -Jx1d -R-170/-130/50/70 -N1 -Df -M -A0/0/1 >> alaska.xy
pscoast -Jx1d -R-170/-130/50/70 -N2 -Df -M -A0/0/1 >> alaska.xy

cat << END > alaska_rivers.xy
> GMT_LONLAT
END
pscoast -Jx1d -R-170/-130/50/70 -Ia -Df -M -A0/0/1 >> alaska_rivers.xy


cat << END > alaska_coasts.xy
> GMT_LONLAT
END
pscoast -Jx1d -R-170/-130/50/70 -W3 -Df -M -A0/0/1 >> alaska_coasts.xy

