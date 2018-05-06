! MODULE: SOLVER
!> @author
!> Kara Abdelhadi, Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet de résoudre un problème d'éléments fini.
MODULE SOLVER
USE SKYLINE
IMPLICIT NONE


CONTAINS
	!> Fonction A(x) du problème donné
	!! \param x Réel
	!! \return A(x)
	FUNCTION A(x)
		REAL :: x
		
		REAL :: A
		
		A = SIN(2.*x) * x
	END FUNCTION A
	
	!> Fonction F(x) du problème donné
	!! \param x Réel
	!! \return F(x)
	FUNCTION F(x)
		REAL :: x
		
		REAL :: F
		
		F = COS(3.*x)
		
	END FUNCTION F
	
	!> Calcule si une valeur est dans un intervalle donné
	!! \param x Valeur
	!! \param xmin Borne inf de l'intervalle
	!! \param xmax Borne max de l'intervalle
	!! \return LOGICAL
	FUNCTION IS_IN_INTERVAL(x, xmin, xmax) RESULT(B)
		REAL :: x
		REAL :: xmin
		REAL :: xmax
		
		LOGICAL :: B
		
		B = .FALSE.
		
		IF ( x<=xmax .AND. x >= xmin ) THEN
			B = .TRUE.
		ENDIF
	END FUNCTION IS_IN_INTERVAL
	
	!> Fonction de base
	!! \param i Indice de la fonction
	!! \param x Valeur d'évaluation
	!! \param h Pas de discrétisation
	!! \param Xdiscret Maillage
	!! \return Evaluation de \f$PHI_i\f$ en x
	FUNCTION PHI(i, x, h, Xdiscret)
		INTEGER :: i
		REAL :: x
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		
		REAL :: PHI
		
		INTEGER :: N
		
		N = SIZE(Xdiscret)
		PHI = 0.
		
		IF (i .EQ. N) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(N-1), Xdiscret(N)) ) THEN
				PHI = (Xdiscret(N) - x) / h
			ENDIF
		ENDIF
		
		IF (i .EQ. 1) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(N-1), Xdiscret(N)) ) THEN
				PHI = (Xdiscret(2) - x) / h
			ENDIF
		ENDIF
		
		IF (i>1 .AND. i<N) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(i - 1), Xdiscret(i))) THEN
				PHI = (x - Xdiscret(i - 1)) / h
			ELSE IF (IS_IN_INTERVAL(x, Xdiscret(i), Xdiscret(i+1))) THEN
				PHI = (Xdiscret(i + 1) - x) / h
			ENDIF
		ENDIF
	END FUNCTION
	
	!> Dérivée de la fonction de base
	!! \param i Indice de la fonction
	!! \param x Valeur d'évaluation
	!! \param h Pas de discrétisation
	!! \param Xdiscret Maillage
	!! \return Evaluation de la dérivée de \f$PHI_i\f$ en x
	FUNCTION PHI_DER(i, x, h, Xdiscret)
		INTEGER :: i
		REAL :: x
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		
		REAL :: PHI_DER
		
		INTEGER :: N
		
		N = SIZE(Xdiscret)
		PHI_DER = 0.
		
		IF (i .EQ. N) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(N-1), Xdiscret(N)) ) THEN
				PHI_DER = -1. / h
			ENDIF
		ENDIF
		
		IF (i .EQ. 1) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(N-1), Xdiscret(N)) ) THEN
				PHI_DER = -1. / h
			ENDIF
		ENDIF
		
		IF (i>1 .AND. i<N) THEN
			IF (IS_IN_INTERVAL(x, Xdiscret(i - 1), Xdiscret(i))) THEN
				PHI_DER = 1. / h
			ELSE IF (IS_IN_INTERVAL(x, Xdiscret(i), Xdiscret(i+1))) THEN
				PHI_DER = -1. / h
			ENDIF
		ENDIF
	END FUNCTION PHI_DER
	
	!> Calcul de l'intégrale du terme bilinéaire
	!! \param binf Borne inf
	!! \param bsup Borne sup
	!! \param N Précision de l'intégrale
	!! \return Valeur de l'intégrale
	FUNCTION trapzA(binf,bsup,N) RESULT(y)
		integer :: i, N
		real ::bsup,binf,h,y
		
		y = 0.
	
		h = (bsup-binf)/N
	
		do i=1,N
			y=y+A(binf + i*h)	
		enddo
	
		y=y+(A(binf)+A(bsup))*0.5
		y=y*h
	end function trapzA
	
	!> Calcul de l'intégrale du terme Linéaire
	!! \param binf Borne inf
	!! \param bsup Borne sup
	!! \param N Précision de l'intégrale
	!! \param ix Indice de la fonction PHI
	!! \param hx Pas du maillage
	!! \return Valeur de l'intégrale
	FUNCTION trapzF(binf,bsup,N, ix, hx, Xdiscret) RESULT(y)
		integer :: i, N
		real ::bsup,binf,h,y
		integer :: ix
		REAL :: hx
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		
		REAL :: x
		y = 0.
	
		h = (bsup-binf)/N
	
		do i=1,N
			x = binf + i*h
			y=y+F(x)*PHI(ix, x, hx, Xdiscret)
		enddo
	
		y=y+(F(binf)*PHI(ix, binf, hx, Xdiscret)+F(bsup)*PHI(ix, bsup, hx, Xdiscret))*0.5
		y=y*h
	end function trapzF
	
	!> Evaluation du terme bilinéaire
	!! \param i Indice Ligne
	!! \param j Indice colonne
	!! \param h Pas du maillage
	!! \param Xdiscret Maillage
	!! \param Nitg Précision intégration
	!! \return 
	FUNCTION B(i, j, h, Xdiscret, Nitg)
		INTEGER :: i
		INTEGER :: j
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		INTEGER :: Nitg
		
		REAL :: B
		
		INTEGER :: N
		
		! PHI(i, x, h, Xdiscret)
		N = SIZE(Xdiscret)
		
		IF ( i .EQ. j) THEN
			IF ( i .EQ. N) THEN
				B = 1./ ( h**2) * trapzA(Xdiscret(N-1), Xdiscret(N), Nitg)
			ELSE IF ( i .EQ. 1) THEN
				B = 1./ ( h**2) * trapzA(Xdiscret(1), Xdiscret(2), Nitg)
			ELSE
				B = 1./ ( h**2) * trapzA(Xdiscret(i - 1), Xdiscret(i + 1), Nitg)
			ENDIF
		ELSE IF ( j .EQ. i+1 ) THEN
			B = -1./(h**2) * trapzA(Xdiscret(i), Xdiscret(i+1), Nitg)
		ELSE IF ( j .EQ. i-1 ) THEN
			B = -1./(h**2) * trapzA(Xdiscret(i-1), Xdiscret(i), Nitg)
		ENDIF
	
	END FUNCTION B
	
	!> Evaluation du terme Linéaire
	!! \param i Indice ligne
	!! \param h Pas du maillage
	!! \param Xdiscret Maillage
	!! \param Nitg Précision intégration
	!! \param alpha1 Parametre du problème
	!! \return 
	FUNCTION L(i, h, Xdiscret, Nitg, alpha1)
		INTEGER :: i
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		INTEGER :: Nitg
		REAL :: alpha1
		
		REAL :: L
		
		INTEGER :: N
		
		! PHI(i, x, h, Xdiscret)
		N = SIZE(Xdiscret)
		
		IF ( i .EQ. 1 ) THEN
			L = trapzF(Xdiscret(1),Xdiscret(2),Nitg, i, h, Xdiscret)
		ELSE IF (i .EQ. N ) THEN
			L = trapzF(Xdiscret(N-1),Xdiscret(N),Nitg, i, h, Xdiscret) + alpha1*A(1.)*PHI(N, 1., h, Xdiscret)
		ELSE
			L = trapzF(Xdiscret(i-1),Xdiscret(i+1),Nitg, i, h, Xdiscret)
		ENDIF
	END FUNCTION L
	
	!> Calcul la matrice de rigidité
	!! \param h Pas du maillage
	!! \param Xdiscret Maillage
	!! \param Nitg Précision intégration
	!! \return Matrice SKYLINE 
	FUNCTION CREATE_STIFFNESS_MATRIX(h, Xdiscret, Nitg) RESULT(K)
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		INTEGER :: Nitg
		
		TYPE(SKYLINE_MATRIX) :: K
		
		INTEGER :: N
		INTEGER :: elem
		INTEGER :: i
		
		N = SIZE(Xdiscret)
		
		K = CREATE_SKYLINE_MATRIX(N, N)
		
		elem = 4
		IF (N>2) THEN
			elem = elem + (N-2)*3
		ENDIF
		
		DEALLOCATE(K%m_S)
		ALLOCATE(K%m_S(elem))
		
		K%m_S(1) = B(1, 1, h, Xdiscret, Nitg)
		K%m_S(2) = B(1, 2, h, Xdiscret, Nitg)
		K%m_Pmoins(1) = 1
		K%m_Pplus(1) = 2
		
		elem = 3
		DO i = 2, N-1, 1
		
			K%m_S(elem)     = B(i, i - 1, h, Xdiscret, Nitg)
			K%m_S(elem + 1) = B(i, i, h, Xdiscret, Nitg)
			K%m_S(elem + 2) = B(i, i + 1, h, Xdiscret, Nitg)
			K%m_Pmoins(i) = i - 1
			K%m_Pplus(i) = i + 1
		
			elem = elem + 3
		ENDDO
		
		K%m_Pmoins(N) = N-1
		K%m_Pplus(N) = N
		K%m_S(elem)     = B(N, N-1, h, Xdiscret, Nitg)
		K%m_S(elem + 1) = B(N, N, h, Xdiscret, Nitg)
	END FUNCTION CREATE_STIFFNESS_MATRIX
	
	!> Calcule du second membre
	!! \param h Pas du maillage
	!! \param Xdiscret Maillage
	!! \param Nitg Précision intégration
	!! \param alpha1 Parametre du problème
	!! \return Matrice SKYLINE 
	FUNCTION CREATE_F_MATRIX(h, Xdiscret, Nitg, alpha1) RESULT(F)
		REAL :: h
		REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		INTEGER :: Nitg
		REAL :: alpha1
		
		REAL , DIMENSION(:), ALLOCATABLE :: F
		
		INTEGER :: N
		INTEGER :: i
		
		N = SIZE(Xdiscret)
		ALLOCATE(F(N))
		
		DO i = 1, N, 1
			F(i) = L(i, h, Xdiscret, Nitg, alpha1)
		ENDDO
	END FUNCTION CREATE_F_MATRIX
	
	!> Résouds le système triangulaire inférieur
	!! \param LT Matrice SKYLINE
	!! \param BB second membre
	!! \return Vecteur Solution
	FUNCTION SOLVE_TRIG_DOWN(LT, BB) RESULT(Y)
		TYPE(SKYLINE_MATRIX) :: LT
		REAL , DIMENSION(:), ALLOCATABLE :: BB

		REAL , DIMENSION(:), ALLOCATABLE ::  Y
		
		INTEGER :: N, i, j
		REAL :: tmp
		
		N = SIZE(BB)
		
		ALLOCATE(Y(N))
		
		DO i = 1, N, 1
			tmp = 0.
			DO j = 1, i - 1, 1
				tmp = tmp + LT%GET(j,i) *Y(j) 
			ENDDO
			Y(i) = (BB(i) - tmp) / LT%GET(i,i)
		ENDDO
	
	END FUNCTION SOLVE_TRIG_DOWN
	
	!> Résouds le système triangulaire supérieur
	!! \param LT Matrice SKYLINE
	!! \param BB second membre
	!! \return Vecteur Solution 
	FUNCTION SOLVE_TRIG_UP(LT, BB) RESULT(Y)
		TYPE(SKYLINE_MATRIX) :: LT
		REAL , DIMENSION(:), ALLOCATABLE :: BB

		REAL , DIMENSION(:), ALLOCATABLE ::  Y
		
		INTEGER :: N, i, j
		REAL :: tmp
		
		N = SIZE(BB)
		
		ALLOCATE(Y(N))
		
		DO i = N, 1, 1
			tmp = 0.
			DO j = i+1, N, 1
				tmp = tmp + LT%GET(j,i) *Y(j) 
			ENDDO
			Y(i) = (BB(i) - tmp) / LT%GET(i,i)
		ENDDO
	
	END FUNCTION SOLVE_TRIG_UP
	
	!> Résouds le problème \f$ LT*Y = BB\f$
	!! \param LT Matrice SKYLINE
	!! \param BB second membre
	!! \return Vecteur Solution 
	FUNCTION SOLVE(LT, BB) RESULT(Y)
		TYPE(SKYLINE_MATRIX) :: LT
		REAL , DIMENSION(:), ALLOCATABLE :: BB

		REAL , DIMENSION(:), ALLOCATABLE ::  Y
		
		Y = SOLVE_TRIG_DOWN(LT, BB)
		Y = SOLVE_TRIG_UP(LT, Y)
	END FUNCTION SOLVE
	
	!> Enregistre un vecteur solution dans un fichier
	!! \param Xdiscret Discrétisation espace
	!! \param U Vecteur solution
	!! \param filename Nom du fichier
	SUBROUTINE EXPORT_SOLUTION(Xdiscret, U, filename)
	        REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret
		REAL , DIMENSION(:), ALLOCATABLE :: U
		character(len=*), intent(in)    :: filename
		
		INTEGER, PARAMETER      :: out_unit = 20
                INTEGER                 :: i
                
                
                open (unit=out_unit,file=filename,action="write",status="replace")
                
                WRITE(out_unit,*) "x,y"
                DO i =1, SIZE(U),1
                        WRITE(out_unit,*) Xdiscret(i),",", U(i)
                ENDDO
                CLOSE (out_unit)
	END SUBROUTINE EXPORT_SOLUTION
END MODULE SOLVER
