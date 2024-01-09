      subroutine eb2mp(mk,n1,n0,sn0,diff,btmn,btva,btdist,ord,grid,btp1)
      integer n1,n0,ord(2,n1+n0),grid(2,n1)
      double precision mk(2,n1+n0),sn0,diff,btmn,btva,btdist(-n0:n0),
     &     btp1(-(n1-1):(n1-1),2)

      external mnprob
      double precision mnprob
      
      integer r,i,j,k,d,drg(2),tmp,a,lg(2),mink,a1,a2
      double precision prb,tri(3)
      logical lgcl,lgcl1

      r=n1-ceiling(dble(n1)*sn0)+1
      diff=0.0d0
      do k=1,2
         do i=1,n1+n0
            ord(k,i)=i
         enddo
         drg(1)=2
         drg(2)=n1
         do d=1,2
            if(d .eq. 2) then
               drg(1)=n1+2
               drg(2)=n1+n0
            endif
            do i=drg(1),drg(2)
               j=i
               lgcl=j .gt. drg(1)-1
               if(lgcl) lgcl=mk(k,ord(k,j-1)) .gt. mk(k,ord(k,j))
               do while(lgcl)
                  tmp=ord(k,j)
                  ord(k,j)=ord(k,j-1)
                  ord(k,j-1)=tmp
                  j=j-1
                  lgcl=j .gt. drg(1)-1
                  if(lgcl) lgcl=mk(k,ord(k,j-1)) .gt. mk(k,ord(k,j))
               enddo
            enddo
         enddo
         tri(1)=mk(k,ord(k,r))
         i=n1+1
         lgcl=mk(k,ord(k,i)) .lt. tri(1)
         do while(lgcl)
            diff=diff+dble(3-2*k)
            i=i+1
            lgcl=i .le. n1+n0
            if(lgcl) lgcl=mk(k,ord(k,i)) .lt. tri(1)
         enddo
         lg(k)=0
         i=0
         j=n1+1
         do while(i .lt. n1)
            i=i+1
            lgcl=j .le. n1+n0
            if(lgcl) lgcl=mk(k,ord(k,j)) .lt. mk(k,ord(k,i))
            do while(lgcl)
               j=j+1
               lgcl=j .le. n1+n0
               if(lgcl) lgcl=mk(k,ord(k,j)) .lt. mk(k,ord(k,i))
            enddo
            if(j .gt. n1+n1) then
               i=n1
            else
               i=i+1
               lgcl=i .le. n1
               if(lgcl) lgcl=mk(k,ord(k,i)) .le. mk(k,ord(k,j))
               do while(lgcl)
                  i=i+1
                  lgcl=i .le. n1
                  if(lgcl) lgcl=mk(k,ord(k,i)) .le. mk(k,ord(k,j))
               enddo
               i=i-1
            endif
            lg(k)=lg(k)+1
            grid(k,lg(k))=i
         enddo
      enddo
      diff=diff/dble(n0)

      a=0
      a1=0
      a2=0
      k=1
      do i=1,3
         tri(i)=0.0d0
      enddo
      do i=-n0,n0
         btdist(i)=0.0d0
      enddo
      d=1
      mink=min(lg(1),lg(2))
      lgcl=k .eq. 1
      if(.not. lgcl) then
         lgcl=k .le. mink
         if(lgcl) lgcl=btp1(0,3-d) .lt. 1.0d0-1.0d-10
      endif
      do while(lgcl)
         i=1
         if(k .gt. 1) i=grid(1,k-1)+1
         do j=i,grid(1,k)
            lgcl1=k .eq. 1
            if(.not. lgcl1)
     &           lgcl1=mk(2,ord(1,j)) .gt. mk(2,ord(2,grid(2,k-1)))
            if(lgcl1) lgcl1=mk(2,ord(1,j)) .le. mk(2,ord(2,grid(2,k)))
            if(lgcl1) a=a+1
         enddo
         a=a+a1+a2
         btp1(0,d)=mnprob(r,n1,a,grid(1,k)-a,grid(2,k)-a)
         drg(1)=a
         i=k+1
         do while(i .le. lg(1))
            do j=grid(1,i-1)+1,grid(1,i)
               if(mk(2,ord(1,j)) .le. mk(2,ord(2,grid(2,k))))
     &              drg(1)=drg(1)+1
            enddo
            if(i .eq. k+1) a1=drg(1)-a
            btp1(i-k,d)=mnprob(r,n1,drg(1),grid(2,k)-drg(1),
     &           grid(1,i)-drg(1))
            if(drg(1) .eq. grid(2,k)) then
               do j=i+1,lg(1)
                  btp1(j-k,d)=btp1(i-k,d)
               enddo
               i=lg(1)
            endif
            i=i+1
         enddo
         drg(2)=a
         i=k+1
         do while(i .le. lg(2))
            do j=grid(2,i-1)+1,grid(2,i)
               if(mk(1,ord(2,j)) .le. mk(1,ord(1,grid(1,k))))
     &              drg(2)=drg(2)+1
            enddo
            if(i .eq. k+1) a2=drg(2)-a
            btp1(k-i,d)=mnprob(r,n1,drg(2),grid(1,k)-drg(2),
     &           grid(2,i)-drg(2))
            if(drg(2) .eq. grid(1,k)) then
               do j=i+1,lg(2)
                  btp1(k-j,d)=btp1(k-i,d)
               enddo
               i=lg(2)
            endif
            i=i+1
         enddo
         
         prb=btp1(0,d)
         if(k .gt. 1) prb=prb+btp1(0,3-d)-btp1(-1,3-d)-btp1(1,3-d)
         if(prb .gt. 1.0d-10) call dspec(mk,n1,n0,ord,grid(1,k),
     &        grid(2,k),tri,btdist,prb)
         do i=k+1,lg(1)
            prb=btp1(i-k,d)-btp1(i-k-1,d)
            if(k .gt. 1) prb=prb+btp1(i-k,3-d)-btp1(i-k+1,3-d)
            if(prb .gt. 1.0d-10) call dspec(mk,n1,n0,ord,grid(1,i),
     &           grid(2,k),tri,btdist,prb)
         enddo
         do j=k+1,lg(2)
            prb=btp1(k-j,d)-btp1(k-j+1,d)
            if(k .gt. 1) prb=prb+btp1(k-j,3-d)-btp1(k-j-1,3-d)
            if(prb .gt. 1.0d-10) call dspec(mk,n1,n0,ord,grid(1,k),
     &           grid(2,j),tri,btdist,prb)
         enddo
         
         k=k+1
         d=3-d
         lgcl=k .eq. 1
         if(.not. lgcl) then
            lgcl=k .le. mink
            if(lgcl) lgcl=btp1(0,3-d) .lt. 1.0d0-1.0d-10
         endif
      enddo

      btmn=tri(1)
      btva=tri(2)-tri(1)**2.0d0+tri(3)

      end

      subroutine dspec(mk,n1,n0,ord,i,j,tri,btdist,prb)
      integer n1,n0,ord(2,n1+n0),i,j
      double precision mk(2,n1+n0),tri(3),btdist(-n0:n0),prb

      integer a(2),b(2),b1,b2,k,p
      double precision dtmp1,dtmp2,la(3)
      logical lgcl

      b(1)=i
      b(2)=j
      do k=1,2
         a(k)=0
         p=n1+1
         lgcl=p .le. n1+n0
         if(lgcl) lgcl=mk(k,ord(k,p)) .lt. mk(k,ord(k,b(k)))
         do while(lgcl)
            if(mk(3-k,ord(k,p)) .ge. mk(3-k,ord(3-k,b(3-k))))
     &           a(k)=a(k)+1
            p=p+1
            lgcl=p .le. n1+n0
            if(lgcl) lgcl=mk(k,ord(k,p)) .lt. mk(k,ord(k,b(k)))
         enddo
      enddo
      tri(1)=tri(1)+dble(a(1)-a(2))/dble(n0)*prb
      tri(2)=tri(2)+(dble(a(1)-a(2))/dble(n0))**2.0d0*prb
      tri(3)=tri(3)+dble(a(1)*(n0-a(1))+a(2)*(n0-a(2))+2.0d0*a(1)*a(2))
     &     /dble(n0)**3.0d0*prb

      if(a(1) .eq. n0 .or. a(2) .eq. n0 .or. a(1)+a(2) .eq. 0) then
         btdist(a(1)-a(2))=btdist(a(1)-a(2))+prb
         return
      endif

      if(a(1) .eq. 0 .or. a(2) .eq. 0 .or. a(1)+a(2) .eq. n0) then
         if(a(1) .eq. 0) then
            la(1)=dlog(dble(a(2)))
            la(2)=dlog(dble(n0-a(2)))
         else
            if(a(2) .eq. 0) then
               la(1)=dlog(dble(a(1)))
               la(2)=dlog(dble(n0-a(1)))
            else
               la(1)=dlog(dble(a(1)))
               la(2)=dlog(dble(a(2)))
            endif
         endif
         do b1=n0,0,-1
            if(b1 .eq. n0) then
               dtmp1=dble(b1)*la(1)-dble(n0)*dlog(dble(n0))
            else
               dtmp1=dtmp1+dlog(dble(b1+1))-dlog(dble(n0-b1))
     &              -la(1)+la(2)
            endif
            if(a(1) .eq. 0) then
               btdist(-b1)=btdist(-b1)+dexp(dtmp1)*prb
            else
               if(a(2) .eq. 0) then
                  btdist(b1)=btdist(b1)+dexp(dtmp1)*prb
               else
                  btdist(b1*2-n0)=btdist(b1*2-n0)+dexp(dtmp1)*prb
               endif
            endif
         enddo
         return
      endif
      
      la(1)=dlog(dble(a(1)))
      la(2)=dlog(dble(a(2)))
      la(3)=dlog(dble(n0-a(1)-a(2)))
      do b1=n0,0,-1
         if(b1 .eq. n0) then
            dtmp1=dble(b1)*la(1)-dble(n0)*dlog(dble(n0))
         else
            dtmp1=dtmp1+dlog(dble(b1+1))-dlog(dble(n0-b1))
     &           -la(1)+la(3)
         endif
         do b2=0,n0-b1
            if(b2 .eq. 0) then
               dtmp2=dtmp1
            else
               dtmp2=dtmp2-dlog(dble(b2))
     &              +dlog(dble(n0-b1-b2+1))+la(2)-la(3)
            endif
            btdist(b1-b2)=btdist(b1-b2)+dexp(dtmp2)*prb
         enddo
      enddo
      
      end

      double precision function mnprob(r,size,np1,np2,np3)
      integer r,size,np1,np2,np3

      integer i,j,k
      double precision a,b,c,d,e,f

      mnprob=0.0d0
      if(np2 .eq. size .or. np3 .eq. size) return
      
      mnprob=1.0d0
      if(np1 .eq. size) return

      do i=0,r-1
         a=0.0d0
         if(i .eq. 0) then
            b=dlog(1.0d0-dble(np1)/dble(size))*dble(size)
            if(min(np2,np3) .gt. 0) then
               c=0.0d0
               do j=1,r
                  c=c+dlog(dble(size-r+j))-dlog(dble(j))
               enddo
            endif
         else
            b=b+dlog(dble(size-i+1))-dlog(dble(i))
     &           +dlog(dble(np1))-dlog(dble(size-np1))
            if(min(np2,np3) .gt. 0)
     &           c=c-dlog(dble(size-i+1))+dlog(dble(r-i+1))
         endif
         mnprob=mnprob-dexp(b)
         if(min(np2,np3) .gt. 0) then
            do j=r-i,size-r
               if(j .eq. r-i) then
                  d=c+dlog(dble(np2))*dble(j)
     &                 +dlog(dble(size-np1-np2))*dble(size-i-j)
     &                 -dlog(dble(size-np1))*dble(size-i)
               else
                  d=d-dlog(dble(j))+dlog(dble(size-i-j+1))
     &                 +dlog(dble(np2))-dlog(dble(size-np1-np2))
               endif
               e=1.0d0
               if(np1+np2+np3 .lt. size) then
                  f=(dlog(dble(size-np1-np2-np3))
     &                 -dlog(dble(size-np1-np2)))*dble(size-i-j)
                  e=1.0d0-dexp(f)
                  do k=1,r-i-1
                     f=f+dlog(dble(size-i-j-k+1))-dlog(dble(k))
     &                    +dlog(dble(np3))-dlog(dble(size-np1-np2-np3))
                     e=e-dexp(f)
                  enddo
               endif
               a=a+dexp(d)*e
            enddo
         endif
         mnprob=mnprob+dexp(b)*a
      enddo
      
      end
