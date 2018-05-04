! MODULE: CHOLESKY
!> @author
!> Kara Abdelhadi, Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet de calculer la décomposition de Cholesky d'une matrice définie positive.

MODULE CHOLESKY_MOD
USE SKYLINE
IMPLICIT NONE

CONTAINS

	!> Calcule la décompostion de Cholesky d'une matrice de type SKYLINE.
	!! \param M Matrice SKYLINE
	!! \return Matrice SKYLINE triangulaire supérieure contenant la décomposition
	FUNCTION CHOLESKY(M) RESULT(L)
   		TYPE(SKYLINE_MATRIX) :: M
        	TYPE(SKYLINE_MATRIX) :: L
            
            	REAL , DIMENSION(:), ALLOCATABLE :: Stmp
        	
        	REAL :: diag
            	REAL :: tmp
        	
        	INTEGER :: i
        	INTEGER :: j
        	INTEGER :: k
            	INTEGER :: ofs
        	! Allocation mémoire
        	L = CREATE_SKYLINE_MATRIX(M%m_m, M%m_m)
        	
            	! Itére sur chaque colonne
            	DO i = 1, M%m_m, 1
                
                	!Calcule de l'élement diagonale
                	diag = 0.
                	tmp = 0.
                	
                	DO k = 1, i - 1, 1
                	    	tmp = tmp + L%GET(k, i)**2
                	ENDDO
                	diag = SQRT(M%GET(i,i) - tmp)
                	
                	! P+ = n 
                	L%m_Pplus(i) = M%m_n
                	! P- = i pour la ligne i
                	L%m_Pmoins(i) = i
                
                	! Alloue la mémoire nécessaire
                	ALLOCATE( Stmp(SIZE(L%m_S) + L%m_Pplus(i) - L%m_Pmoins(i) + 1))
                	Stmp = 0.
                	Stmp(1:SIZE(L%m_S)) = L%m_S(1:SIZE(L%m_S))
                	ofs = SIZE(L%m_S) + 1
                	CALL MOVE_ALLOC(Stmp, L%m_S)
                	!! WRITE(*,*) "Line ", i, " Data ", ofs
                	! Assigne l'élément diagonale
                	L%m_S(ofs) = diag
                
                	! On utilise la matrice triangulaire supérieure pour stocker L
                	! De cette façon, on accede a la matrice de maniere linéaire
                	DO j = i + 1, M%m_n, 1 
                    		tmp = 0.
                    		DO k = 1, i - 1, 1
                    		    tmp = tmp + L%GET(k, i)*L%GET(k, j)
                    		ENDDO
                    
                    		L%m_S(ofs + j - i) = (M%GET(i, j) - tmp) / diag
                	ENDDO
            ENDDO
	END FUNCTION CHOLESKY

END MODULE CHOLESKY_MOD
