! MODULE: SKYLINE_MATRIX
!> @author
!> Kara Abdelhadi, Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet de stocker une matrice sous la forme d'une matrice Skyline (profil)
MODULE SKYLINE
IMPLICIT NONE


TYPE SKYLINE_MATRIX
	INTEGER :: m_m ! Line
	INTEGER :: m_n ! Column
	
	INTEGER , DIMENSION(:), ALLOCATABLE     :: m_Pmoins
	INTEGER , DIMENSION(:), ALLOCATABLE     :: m_Pplus
	REAL , DIMENSION(:), ALLOCATABLE     	:: m_S

	CONTAINS
		PROCEDURE       ::  GET, DISPLAY, PROD_MAT_VEC, DISPLAY_DEBUG
END TYPE SKYLINE_MATRIX


CONTAINS
	!> Créé une matrice SKYLINE vide de dimension \f$(m,n)\f$
	!! \param m Nombre de ligne
	!! \param n Nombre de colonne
	!! \return Matrice SKYLINE
	FUNCTION CREATE_SKYLINE_MATRIX(m, n)
		INTEGER, intent(in) :: m
        	INTEGER, intent(in) :: n
		
		TYPE(SKYLINE_MATRIX) CREATE_SKYLINE_MATRIX
		
		CREATE_SKYLINE_MATRIX%m_m = m
        	CREATE_SKYLINE_MATRIX%m_n = n
        	        
        	ALLOCATE( CREATE_SKYLINE_MATRIX%m_Pmoins(m) )
        	ALLOCATE( CREATE_SKYLINE_MATRIX%m_Pplus(m) )
        	
        	ALLOCATE( CREATE_SKYLINE_MATRIX%m_S(0) )
	END FUNCTION CREATE_SKYLINE_MATRIX
	
	!> Créé une matrice SKYLINE depuis les données contenus dans un fichier
	!! \param filename Nom du fichier
	!! \return Matrice SKYLINE
	FUNCTION CREATE_SKYLINE_MATRIX_FROM_FILE( filename ) RESULT(M)
		character(len=*), intent(in)    :: filename
		TYPE (SKYLINE_MATRIX) M
		
		INTEGER, PARAMETER      :: out_unit = 20
                INTEGER                 :: tmp
		
		open (unit = out_unit, file = filename)
                
                READ(out_unit,*) M%m_m
                READ(out_unit,*) M%m_n
                
                ALLOCATE( M%m_Pmoins(M%m_m) )
        	ALLOCATE( M%m_Pplus(M%m_m) )
                
                READ(out_unit,*) M%m_Pmoins
                READ(out_unit,*) M%m_Pplus
                
                READ(out_unit,*) tmp
        
                ALLOCATE(M%m_S(tmp))
                READ(out_unit,*) M%m_S
                
                CLOSE(out_unit)
	END FUNCTION CREATE_SKYLINE_MATRIX_FROM_FILE
	
	!> Sauvegarde une matrice SKYLINE dans un fichier.
	!! \param M Matrice à sauvegarder
	!! \param filename Nom du fichier
	SUBROUTINE WRITE_SKYLINE_MATRIX_TO_FILE(M, filename )
		character(len=*), intent(in)    :: filename
		TYPE (SKYLINE_MATRIX) M
		
		INTEGER, PARAMETER      :: out_unit = 20
                INTEGER                 :: i
                
                
                open (unit=out_unit,file=filename,action="write",status="replace")
                
                WRITE(out_unit,*) M%m_m
                WRITE(out_unit,*) M%m_n
                
                WRITE(out_unit,*) (M%m_Pmoins(i), i=1,M%m_m)
                
                WRITE(out_unit,*) (M%m_Pplus(i), i=1,M%m_m)
                
                WRITE(out_unit,*) SIZE(M%m_S)
                
                WRITE(out_unit,*) (M%m_S(i), i=1,SIZE(M%m_S))
                
                CLOSE (out_unit)
	END SUBROUTINE WRITE_SKYLINE_MATRIX_TO_FILE

	!> Supprime les allocations mémoires d'une matrice SKYLINE
	!! \param M Matrice SKYLINE
	SUBROUTINE FREE_SKYLINE_MATRIX(M)
		TYPE(SKYLINE_MATRIX) ::  M
        	        
        	M%m_n   = 0
        	M%m_m   = 0
        	        
        	DEALLOCATE( M%m_Pmoins )
        	DEALLOCATE( M%m_Pplus )
        	DEALLOCATE( M%m_S )
	END SUBROUTINE FREE_SKYLINE_MATRIX

	!> Accède aux éléments d'une matrice SKYLINE
	!! \param M Matrice SKYLINE
	!! \param row Ligne de l'élément
	!! \param col Colonne de l'élément
	!! \return Valeur de M(i,j)
	FUNCTION GET(M, row, col)
                CLASS(SKYLINE_MATRIX), intent(in)       ::  M
                INTEGER, intent(in)                     :: row
                INTEGER, intent(in)                     :: col
                
                REAL  :: GET
                
                
                ! Check if the value is in Pmoins(i) , Pplus(i)
                
                IF ( M%m_Pmoins(row) > col .OR. M%m_Pplus(row) < col ) THEN
                	GET = 0.
                	RETURN
                ENDIF
                
                GET = M%m_S(POSITION_ON_S(M, row, col))
                
        END FUNCTION GET
        
        !> Calcule la position d'un élément dans le tableau stockant les données
	!! \param M Matrice SKYLINE
	!! \param row Ligne de l'élément
	!! \param col Colonne de l'élément
	!! \return Indice de l'élément dans le tableau de donnée
        FUNCTION POSITION_ON_S(M, row, col) RESULT(K)
        	CLASS(SKYLINE_MATRIX), intent(in)       ::  M
                INTEGER, intent(in)                     :: row
                INTEGER, intent(in)                     :: col
        	
        	INTEGER :: K
        	
        	INTEGER :: beg
        	INTEGER :: i
        	
        	IF ( row < 1 .OR. row>M%m_m .or. col < 1 .OR. col>M%m_n) THEN
        		K = 0
        		RETURN
        	ENDIF 
        	
        	! Calcul de le début de la ligne sur S
		beg = 1
		
		DO i=1, row - 1, 1
			beg = beg + NUMBER_ITEM_ON_LINE(M, i)
		ENDDO
		
		K = beg + col - M%m_Pmoins(row) 
        END FUNCTION POSITION_ON_S
        
        !> Calcule le nombre d'élément stocké sur une ligne
	!! \param M Matrice SKYLINE
	!! \param row Indice de ligne
	!! \return  Nombre d'élément stocké sur une ligne
        FUNCTION NUMBER_ITEM_ON_LINE(M, row) RESULT(N)
        	CLASS(SKYLINE_MATRIX), intent(in)       ::  M
                INTEGER, intent(in)                     :: row
                
                INTEGER :: N
                
                N = M%m_Pplus(row) - M%m_Pmoins(row) + 1
        END FUNCTION NUMBER_ITEM_ON_LINE
        
        !> Affiche une matrice SKYLINE
	!! \param M Matrice SKYLINE
        SUBROUTINE DISPLAY(M)
                CLASS(SKYLINE_MATRIX), intent(in) :: M
                
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: tmp_row 
                INTEGER                                 :: i
                INTEGER                                 :: j

                ALLOCATE(tmp_row(M%m_n))
                
                DO i = 1, M%m_m, 1
                        DO j = 1, M%m_n, 1
                                ! Use a buffer to display a line because WRITE return to line after each call 
                                ! Can be fixed using some weird parameters for WRITE
                                tmp_row(j) = GET(M, i, j)
                        ENDDO
                        WRITE(*,*) tmp_row
                ENDDO
                WRITE(*,*) ""
        END SUBROUTINE DISPLAY
        
        !> Affiche des informations de debuguage d'une matrice SKYLINE
	!! \param M Matrice SKYLINE
        SUBROUTINE DISPLAY_DEBUG(M)
                CLASS(SKYLINE_MATRIX), intent(in) :: M
                
                ! Raw access to CSR_MATRIX
                WRITE(*, *) "DEBUG INFO"
                WRITE(*, *) "m_m  ", M%m_m
                WRITE(*, *) "m_n  ", M%m_n
                WRITE(*, *) "P-    :", M%m_Pmoins
                WRITE(*, *) "P+    :", M%m_Pplus
                WRITE(*, *) "S     :", M%m_S
                CALL M%DISPLAY()
        END SUBROUTINE DISPLAY_DEBUG
        
        !> Calcule Le produit entre une matrice et un vecteur
	!! \param M Matrice SKYLINE
	!! \param X Vecteur
	!! \return Vecteur = M*X
        FUNCTION PROD_MAT_VEC(M, X) RESULT(Y)
        	CLASS(SKYLINE_MATRIX) M
        	REAL , DIMENSION(:), ALLOCATABLE :: X
        	
        	REAL , DIMENSION(:), ALLOCATABLE :: Y
        	
        	INTEGER :: i, j, k
        	REAL :: tmp
        	
        	ALLOCATE(Y(M%m_m))
        	k = 1
        	DO i = 1,M%m_m,1
        		tmp = 0.
        		DO j = M%m_Pmoins(i), M%m_Pplus(i),1
        			tmp = tmp + M%m_S(k) * X(j)
        			k= k + 1
        		ENDDO 
        		Y(i) = tmp
        	ENDDO
        END FUNCTION PROD_MAT_VEC
        
END MODULE












