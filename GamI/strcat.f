	character*80 function strcat(n, string, ichoice)
c
c STRCAT concatenates strings of length 80
c
c ichoice = 0 remove trailing blanks and tabs 
c ichoice = 1 remove all blanks and tabs
c
	integer n, ichoice
	character*80 string(n)
	character*1 c
c
	integer i, j, jj, kk, k2, ipass, ii
	character*80 str, str2
	logical flag
c
c Read strings backwards starting from last ...
c
	k2 = 0
	flag = .false.
c
	do 100 i = 1, n
		ii = n + 1 - i
		str = string(ii)
		if (ichoice .eq. 0) ipass = 0
		do 101 j = 1, 80
			jj = 81 - j
c			print *,' ii jj character =',ii,jj,' ',str(jj:jj)
			c = str(jj:jj)
			if ((c .ne. ' ' .and. c .ne. '	') .or.
     &			 (ichoice .eq. 0 .and. ipass .eq. 1)) then
				ipass = 1
				k2 = k2 + 1
				if (k2 .gt. 80) then
					flag = .true.
					k2 = 80
				endif
				str2(k2:k2) = str(jj:jj)
c				print *,' k2 =',k2
				if (k2 .eq. 80) go to 111
			endif
101		continue
100	continue
c
111	continue
c
	if (flag) print *,' STRCAT: over 80 chars'
c
c Invert string
c
	do 200 k = 1, 80
		kk = k2 + 1 - k
		if (kk .gt. 0) then
			strcat(k:k) = str2(kk:kk)
		else
			strcat(k:k) = ' '
		endif
200	continue
c
	return
	end
c
c
c
