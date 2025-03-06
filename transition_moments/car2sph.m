function c = car2sph(l,m,lx,ly,lz)

  Ltotal = lx+ly+lz;
  if (Ltotal ~= l )
      c = 0;
      return;
  end
  
  j = (lx+ly-abs(m))/2 ;
  if (rem (j , 1 ) >0)
      c = 0;
      return;
  end
  
  if (Ltotal == l)
      pref = sqrt(factorial(2*lx)*factorial(2*ly)*factorial(2*lz)*factorial(l)*factorial(l-abs(m)) ... 
          /(factorial(2*l)*factorial(lx)*factorial(ly)*factorial(lz)*factorial(l+abs(m)))) ...
          /(factorial(l)*2^l) ;

      sumval = 0;
      i = 0;
      while (i <= (l-abs(m))/2)
          sumsumval = 0 ;
          for k = 0 : j
              if m>=0
                  ttmmpp = (abs(m)-lx+2*k)/2;
              elseif m<0
                  ttmmpp = -(abs(m)-lx+2*k)/2;
              end

              if ((abs(m)>=(lx-2*k))&((lx-2*k)>=0))
                  absmchooselxm2k = nchoosek(abs(m),lx-2*k);
              else
                  % if (abs(m)<(lx-2*k))
                  absmchooselxm2k = 0;
              end
              absmchooselxm2k;
              
              sumsumval = sumsumval +nchoosek(j,k)*absmchooselxm2k*(-1)^(ttmmpp);
          end
          if ((i < j)|(j<0))
              ichoosej = 0;
          elseif (i>=j)
              ichoosej = nchoosek(i,j);
          end
          sumval = sumval + nchoosek(l,i)*ichoosej*(-1)^i*factorial(2*l-2*i)/factorial(l-abs(m)-2*i)*sumsumval;
          i = i+1;
      end
      c = sumval * pref ; 
      
  end
  return 
      