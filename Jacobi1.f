      program main
      implicit none !Jacobi iterative method, x(k+1) = (I-inv(D))*x(k) +inv(D)*b
      real*8 :: C, w, test, dt, XOLD(400,400), X(400,400), XNEW(400,400)
     $ ,start_time, end_time, elapsed_time
      integer :: I, J, n, m, nsteps, MAX_ITER, step
	PARAMETER (C=10.0D0, w=1.65, n=400, MAX_ITER=1000, nsteps=100)

	
	dt = 4.0D0 / n
	DO I = 1, n
         DO J = 1, n
	      XOLD(I,J) = exp(-2*((-2+dt*I)**2+(-2+dt*J)**2)) !initial condition
		  X(I,J) = 0.0D0 
         END DO
      END DO

	DO step = 1, nsteps
	call CPU_TIME (start_time) !start clock
         X = XOLD
	Do m = 1, MAX_ITER
        DO I= 2, n-1
	     DO J = 2, n-1
           XNEW(I,J) = C/(4*C+1)*(X(I-1,J)+X(I+1,J)+X(I,J-1)+X(I,J+1))
     $     +1/(4*C+1)*XOLD(I,J)		   
	     END DO
         END DO
         IF(SUM(ABS(X-XNEW)).LT.0.000001) THEN
	      exit

         END IF
         X = XNEW

	END DO
	XOLD=X
	print*, ABS(SUM(X-XNEW)),m
      OPEN(01, STATUS='UNKNOWN', FILE='X.dat', ACCESS='APPEND')
             DO I = 1, n
                WRITE(01,46) (X(I,J), J = 1, n)
            END DO
      CLOSE(01)
   46 FORMAT(1X, 220(F10.6,1X))
      call CPU_TIME (end_time)
      elapsed_time = end_time - start_time 
	print*, elapsed_time
	OPEN(02, STATUS='UNKNOWN', FILE='elapsed_time.dat', 
     $	ACCESS='APPEND')
	WRITE(01,*) elapsed_time
	CLOSE(02)
	END DO

	
      end program main