		 subroutine centra(zulu,imin,imax)
         implicit none
         integer imin,imax,i
		 character *80 zulu
		 do 1 i=1,80
		 if (zulu(i:i).ne.' ') goto 11
1		 continue
11		 imin=i
		 do 2 i=1,80
		 if (zulu(81-i:81-i).ne.' ') goto 22
2		 continue
22		 imax=81-i
		 return
		 end
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
