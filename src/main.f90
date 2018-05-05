!> 
!! @mainpage  Utilisation de matrice SKYLINE et de la décomposition de Cholesky pour la résolution d'un problème d'éléments finis
!! \section Présentation
!! On utilise des matrices SKYLINE pour représenter les matrices en mémoires. On peut alors en calculer une décomposition de Cholesky. Enfin, on peut résoudre le problème d'élément finis.
!! \section Cholesky
!! Ici on ne stocke pas la matrice L calculé, mais sa transposée. Cela permet de n'avoir a gérer que des acces en ligne durant les phases d'écriture.
!! \section SOLVER
!! Il faut modifier A(x) et F(x) pour résoudre le problème voulu.
!! \section Compilation
!! Pour compiler le code executer le script *compile.sh*.


PROGRAM MAIN
        USE SKYLINE
        USE SOLVER
        USE CHOLESKY_MOD
        
        INTEGER :: N, i
        TYPE(SKYLINE_MATRIX) :: S_K, S_L
        REAL , DIMENSION(:), ALLOCATABLE :: Xdiscret, S_F, S_U
        
        REAL :: alpha1, h
        INTEGER :: Nitg
        N = 50
        
        
        ! Définition de la discrétisation 
        h = 1./(N - 1)
        ALLOCATE(Xdiscret(N))
        DO i = 1, N, 1
        	Xdiscret(i) = (i - 1.) * h
        ENDDO
        
        alpha1 = 1.
        Nitg = 10
        
        ! Création des matrices décrivant le problème variationnel discret
        S_K = CREATE_STIFFNESS_MATRIX(h, Xdiscret, Nitg)
        S_F = CREATE_F_MATRIX(h, Xdiscret, Nitg, alpha1)
        
        ! Calcul de la décomposition de cholesky de la matrice S_K
        S_L = CHOLESKY(S_K)
        
        ! Résolution du probleme variationnel discret
        S_U = SOLVE(S_K, S_F)
        
        ! Libération mémoire
	DEALLOCATE(S_U)
	DEALLOCATE(S_F)
	DEALLOCATE(Xdiscret)
	
	CALL FREE_SKYLINE_MATRIX(S_K)
	CALL FREE_SKYLINE_MATRIX(S_L)
        
END PROGRAM MAIN





