//Questao 1 - P2 CV916 

// Coordenadas dos nos
Int XY = [	0	  0
		1/12	  0
		2/12     0
		3/12     0
		4/12     0
	 	5/12     0
		6/12     0
		7/12     0
		8/12     0
		9/12     0
		10/12     0
		11/12     0
		1       0]

// Criacao dos elementos
Int Beam_EL = [	1	2	
			2     3
			3	4
			4	5
			5	6
			6	7
			7	8
			8	9
			9	10
			10	11
			11	12
			12	13]	
			  
// Definicao dos graus de liberdade para cada no
Int EQ = [	0 0 1
		0 2 3
		0 4 5
		0 6 7
		0 8 9
		0 10 11
  		0 12 13
		0 14 15
		0 16 17
		0 18 19
		0 20 21
		0 22 23
		0 0 24]

// Locacao dos elementos na matriz
Int Beam_LM = Beam_Make_LM(Beam_EL, EQ) 

// Propriedades das matrizes [A As I E Nu Dn]
Dbl P1 = [1.0 0.0 1.0 1.0 0.0 0.0]
Dbl Prop = [1.0 0.0 1.0 1.0 0.0 0.0] 

// Definicao da matriz de rigidez local
Dbl KE1 = Beam_Make_K(1,XY,Beam_EL,P1,"Global")
Dbl KE2 = Beam_Make_K(2,XY,Beam_EL,P1,"Global") 
Dbl KE3 = Beam_Make_K(3,XY,Beam_EL,P1,"Global")  
Dbl KE4 = Beam_Make_K(4,XY,Beam_EL,P1,"Global")
Dbl KE5 = Beam_Make_K(5,XY,Beam_EL,P1,"Global") 
Dbl KE6 = Beam_Make_K(6,XY,Beam_EL,P1,"Global") 
Dbl KE7 = Beam_Make_K(7,XY,Beam_EL,P1,"Global")
Dbl KE8 = Beam_Make_K(8,XY,Beam_EL,P1,"Global") 
Dbl KE9 = Beam_Make_K(9,XY,Beam_EL,P1,"Global") 
Dbl KE10 = Beam_Make_K(10,XY,Beam_EL,P1,"Global")
Dbl KE11 = Beam_Make_K(11,XY,Beam_EL,P1,"Global") 
Dbl KE12 = Beam_Make_K(12,XY,Beam_EL,P1,"Global") 

// Montagem da matriz de rigidez global
Dbl KE = Zero (24) 
Assemble (KE, KE1, Beam_LM, 1)
Assemble (KE, KE2, Beam_LM, 2)
Assemble (KE, KE3, Beam_LM, 3)
Assemble (KE, KE4, Beam_LM, 4)
Assemble (KE, KE5, Beam_LM, 5)
Assemble (KE, KE6, Beam_LM, 6)
Assemble (KE, KE7, Beam_LM, 7)
Assemble (KE, KE8, Beam_LM, 8)
Assemble (KE, KE9, Beam_LM, 9)
Assemble (KE, KE10, Beam_LM, 10)
Assemble (KE, KE11, Beam_LM, 11)
Assemble (KE, KE12, Beam_LM, 12)

// definicao da matriz de rigidez geométrica local
//KG = Beam_Make_KG(N ,e , XY, Beam_EL, Prop, Ref, Type)
Dbl KC11 = Beam_Make_KG(1,1,XY,Beam_EL,Prop,"Global","C")
Dbl KC22 = Beam_Make_KG(1,2,XY,Beam_EL,Prop,"Global","C")
Dbl KC33 = Beam_Make_KG(1,3,XY,Beam_EL,Prop,"Global","C")
Dbl KC44 = Beam_Make_KG(1,4,XY,Beam_EL,Prop,"Global","C")
Dbl KC55 = Beam_Make_KG(1,5,XY,Beam_EL,Prop,"Global","C")
Dbl KC66 = Beam_Make_KG(1,6,XY,Beam_EL,Prop,"Global","C")
Dbl KC77 = Beam_Make_KG(1,7,XY,Beam_EL,Prop,"Global","C")
Dbl KC88 = Beam_Make_KG(1,8,XY,Beam_EL,Prop,"Global","C")
Dbl KC99 = Beam_Make_KG(1,9,XY,Beam_EL,Prop,"Global","C")
Dbl KC1010 = Beam_Make_KG(1,10,XY,Beam_EL,Prop,"Global","C")
Dbl KC1111 = Beam_Make_KG(1,11,XY,Beam_EL,Prop,"Global","C")
Dbl KC1212 = Beam_Make_KG(1,12,XY,Beam_EL,Prop,"Global","C")

// montagem da matriz de rigidez geométrica global
Dbl KG = Zero (24)
Assemble (KG, KC11, Beam_LM, 1)
Assemble (KG, KC22, Beam_LM, 2)
Assemble (KG, KC33, Beam_LM, 3)
Assemble (KG, KC44, Beam_LM, 4)
Assemble (KG, KC55, Beam_LM, 5)
Assemble (KG, KC66, Beam_LM, 6)
Assemble (KG, KC77, Beam_LM, 7)
Assemble (KG, KC88, Beam_LM, 8)
Assemble (KG, KC99, Beam_LM, 9)
Assemble (KG, KC1010, Beam_LM, 10)
Assemble (KG, KC1111, Beam_LM, 11)
Assemble (KG, KC1212, Beam_LM, 12)

Dbl {Modos,lambda} = Jacobi(KG,KE)

End
