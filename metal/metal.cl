procedure metal (o3wflux, o3wunc, o2flux, o2unc, o3flux, o3unc, hbflux, hbunc, av1, count1)
real	o3wflux			{prompt="Line flux of O3 4363 line"}
real	o3wunc			{prompt="The uncertainity in the O3 4363 line"}
real	o2flux			{prompt="The line flux of O2 3728 line"}
real	o2unc			{prompt="The uncertainity in the O2 3728 line"}
real	o3flux			{prompt="Line flux of O3 4959 line"}
real	o3unc			{prompt="The uncertainity in the O3 4959 line"}
real	hbflux			{prompt="The line flux of Hb 4861 line"}
real	hbunc			{prompt="The uncertainity in the hb 4861 line"}
real	av1			{prompt="The av value from starlight for reddening"}
real	count1			{prompt="The maximum count value"}
# list of input variables, all relevant fluxes and the reddening value

struct	*flist



begin
	real	o3wf, o3wu, o2f, o2u, o3f, o3u, hbf, hbu, v1, v2, v3, v4, n1, n2, n3, n4, num 	#flux variables, then rnd nos
	real	o3wfm, o2fm, o3fm, hbfm								#changing flux from each monte sweep
	real	r1, r2										# radius of box muller
	real	o3temp, o2temp, o3abund, o2abund, metal						#calculated values
	real	wvx, wvy, ax, bx, ao2, ahb, av, ratio						# reddening values
	real	summet, summet2, sum3t, sum3a, sum2a 						# average values
	string tmprnd

	real count, loop, max									# interation values
	
	max = count1
	o3wf=o3wflux
	o3wu=o3wunc
	o2f=o2flux
	o2u=o2unc
	o3f=4*o3flux
	o3u=4*o3unc
	hbf=hbflux	
	hbu=hbunc
	av = av1
	summet = 0
	summet2 = 0
	sum3t = 0
	sum3a = 0
	sum2a = 0
	tmprnd = "tmprnd.txt"
	count = 0

# find reddening ratio using CCM formulae
		wvx = 1 / 0.3728
	      	wvx = (wvx - 1.82)
      		ax = 1.0 + (0.17699 - 0.50447*wvx)*wvx
		ax -= 0.02427*wvx**3
		ax += 0.72085*wvx**4
		ax += 0.01979*wvx**5
		ax -= 0.77530*wvx**6
		ax += 0.32999*wvx**7

		bx = 1.41338*wvx
		bx += 2.28305*wvx**2
		bx += 1.07233*wvx**3
		bx -= 5.38434*wvx**4
		bx -= 0.62251*wvx**5
		bx += 5.30260*wvx**6
		bx -= 2.09002*wvx**7

 		ao2 = (ax + bx/3.1) * av

		wvx = 1 / 0.4861
	      	wvx = (wvx - 1.82)
      		ax = 1.0 + (0.17699 - 0.50447*wvx)*wvx
		ax -= 0.02427*wvx**3
		ax += 0.72085*wvx**4
		ax += 0.01979*wvx**5
		ax -= 0.77530*wvx**6
		ax += 0.32999*wvx**7

		bx = 1.41338*wvx
		bx += 2.28305*wvx**2
		bx += 1.07233*wvx**3
		bx -= 5.38434*wvx**4
		bx -= 0.62251*wvx**5
		bx += 5.30260*wvx**6
		bx -= 2.09002*wvx**7

 		ahb = (ax + bx/3.1) * av	

		ratio = 10 ** (0.4 * (ao2-ahb))
		print(ratio)

		# Create a lot of rnd no. might need to increase no of rnd no for max = 1 case
		delete(tmprnd,verify-)
		urand(seed = indef, nlines = 2*max, ncols = 4, >>tmprnd)

		flist = tmprnd


	while (count != max){
		count += 1
		metal = 0
		o3temp = 0
		o2temp = 0
		o3abund = 0
		o2abund = 0
		#values for each metal calc

		#make rnd no. using urand and box muller
		loop = 0
		while(loop != 1){

			#Moved urand to before loop 22/12
#			delete(tmprnd,verify-)
#			urand(seed = indef, nlines = 1, ncols = 4, >>tmprnd)

#			flist = tmprnd
		
			if (fscan(flist,v1,v2, v3, v4) == 4) {
#				flist = " "
			} else {
				print ("Error w/ random no")
				bye
				}
			v1 = v1 * 2 -1
			v2 = v2 * 2 -1
			v3 = v3 * 2 -1
			v4 = v4 * 2 -1
			r1 = v1**2 + v2**2
			r2 = v3**2 + v4**2
			
			if (r1 < 1 && r1 != 0 && r2 < 1 && r2 != 0) {
				loop = 1
				n1 = v1 * sqrt(-2 * log(r1) / r1)
				n2 = v2 * sqrt(-2 * log(r1) / r1)	

				n3 = v3 * sqrt(-2 * log(r2) / r2)
				n4 = v4 * sqrt(-2 * log(r2) / r2)
				if ((o3wf + n1*(o3wu))<0){
					loop = 0 }
				}
			}
	#print(count)	

	#find o3 abund
		# randomised fluxes with standard deviation
		o3wfm = o3wf + n1*(o3wu)
		o2fm = o2f + n2*(o2u)
		o3fm = o3f + n3*(o3u)
		hbfm = hbf + n4*(hbu)
		
		#print(o3f, o3wf)
		
		# use temden to get temp		
		temden("temperature", flxratio= o3fm/o3wfm, assume = 100, spec = 3, atom = "oxygen",>> "log.txt")
		o3temp = temden.result
		sum3t += o3temp
		#print(o3temp,o3wfm, o3wf, n1, o3wu)
		o2temp = 0.7*o3temp + 3000
		
		#use ionic to get abund
		ionic(atom = "oxygen", spectrum = 3, temperature = o3temp, density = 100.0, wave = 4983, wv_tol = 30, flx = o3fm*100/hbfm, >> "log.txt")

		o3abund = ionic.result
		sum3a += o3abund
#	print (o3abund)

	# find o2 abundance

		ionic(atom = "oxygen", spectrum = 2, temperature = o2temp, density = 100.0, wave = 3728, wv_tol = 3, flx = 1.0464*ratio*o2fm*100/hbfm, >> "log.txt")

		o2abund = ionic.result
		sum2a += o2abund

		metal = 12+log10(o3abund + o2abund)
		summet += metal
		summet2 += metal ** 2
		#print(o3fm, o3f, n3, o3u, >> "wv.txt")			# debug prints for more information
		#print(o3wfm, o2fm, o3fm, hbfm, metal, >> "metal.txt") 	# debug prints for more information


	}
	# average values
	summet = summet / count
	summet2 = sqrt(summet2/ count - (summet)**2)
	sum3t = sum3t / count
	sum3a = sum3a / count
	sum2a = sum2a / count

	print("The average metallicity is; ", summet)
	print("The standard deviation is; ", summet2)

	#debug values
	print("Average o3 temp = ", sum3t)
	print("Average o3 abun = ", sum3a)
	print("Average o2 abun = ", sum2a)
	
	

end
