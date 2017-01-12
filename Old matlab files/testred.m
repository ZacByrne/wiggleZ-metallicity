wv = 0.4861;
wvx = 1 / wv;
wvx = (wvx - 1.82);
ax = 1.0 + (0.17699 - 0.50447*wvx)*wvx;
		ax = ax - 0.02427*wvx^3;
		ax = ax + 0.72085*wvx^4;
		ax = ax + 0.01979*wvx^5;
		ax = ax - 0.77530*wvx^6 ;
		ax = ax + 0.32999*wvx^7;

		bx = 1.41338*wvx;
		bx = bx + 2.28305*wvx^2;
		bx = bx + 1.07233*wvx^3;
		bx = bx - 5.38434*wvx^4;
		bx = bx - 0.62251*wvx^5;
		bx = bx + 5.30260*wvx^6;
		bx = bx - 2.09002*wvx^7;

        
        ax
    
        bx
        wv
        k = ax*3.1 + bx