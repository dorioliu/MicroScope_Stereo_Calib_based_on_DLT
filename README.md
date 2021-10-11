# MicroScope_Stereo_Calib_based_on_DLT
 /*********  
 	* This is an implementation of DLT based camera calibration with Eigen  
 	*  
 	* General,we describe a imaging system of pinhole model as:  
    
 	*             sp = K[R|t]P   (1)   
    
  * s is depth scale coefficient, lower case p     
 	* represents 2d image point (u,v,1), upper P is a space 3d  
 	* point (X,Y,Z,1)  in world coordinate.    
 	*  
 	* usually, in rigid transformation, we descibe a rigid motion  
 	* from P1(X1,Y1,Z1) to P2(X2,Y2,Z2) with rotation transformation  
 	* R and the new camera position is C as P2 = R(P1 - C), also  
 	* sometimes, we hide description of C and introduces a tanslation  
 	* vecton t = -RC ,in this case, rigid transformation can also  
 	* be written down as P2 = RP1 + t.  
 	*  
 	* So:  
 	*     sp = K[R|t]P = [KR, -KRC]P = HP  (2)  
 	*  
 	* which C reperents the camera optical center position in world  
 	* coordinate, H is a preojection matrix with a shape of 3 by 4.  
 	*  
 	* The DLT method takes a set of at least 6 pairs of 3d world  
 	* points and those corrsponding projected 2d iamge points as  
 	* inputs and outputs the K,R parameters what we want to get  
 	* by DLT based calibration. please note that all 3d world points  
 	* must not space on a space planar.  
 	*  
 	* So, how to evaluate H and how to do decompsition about H  
 	* to get K and R ?  
 	* in general, we calculate H by build over-determined equation  
 	* based on the relation between p and P.  

  
 	* [u]    [ h11, h12, h13, h14 ]  
 	* [v] =  [ h21, h22, h23, h24 ] * [X, Y,Z,1]^T   (3)      
 	* [1]    [ h31, h32, h33, h34 ]  
 	 
 	*  
 	* thus, we get two quations:  
 	*  
 	* u = (hllX + h12Y + h13Z + 1)/(h31X + h31Y + h33Z + 1)      (4)     
 	* v = (h2lX + h22Y + h23Z + 1)/(h31X + h31Y + h33Z + 1)  
 	*  
 	* now we have I pairs of that point (ui,vi,1) and (Xi,Yi,Zi,1)  
 	* each point pair provide two equation.  
 	*  
 	* based on (4) we get a linear algebraic equation:  
 	*  
 	* [ -X1 -Y1 -Z1 -1 0 0 0 0 x1X1 x1Y1 x1Z1 x1 ]  
 	* [ 0 0 0 0 -X1 -Y1 -Z1 -1 y1X1 y1Y1 y1Z1 y1 ]  
 	*                      .....  
 	* [ -Xi -Yi -Zi -1 0 0 0 0 xiXi xiY1 xiZi xi ]    * [h11, h12, h13, h14,h21, h22, h23, h24,h31, h32, h33, h34 ]^T = 0 (5)  
 	* [ 0 0 0 0 -Xi -Yi -Zi -1 yiXi yiYi yiZi yi ]  
 	*                     .....  
 	* [ -XI -YI -ZI -1 0 0 0 0 xIXI xIYI xIZI xI ]  
 	* [ 0 0 0 0 -XI -YI -ZI -1 yIXI yIYI yIZI yI ]  
 	*  
 	* this linear equation is a traditional Ax = 0 problem.  
 	* the sulution of x can be evaluated by SVD decomposition  
 	* (sigular value decomposition)  
 	* at that case:  
 	* A = USV^T, U: 2Ix12, V: 12x12, S: 12x12 is diagonal matrix with descendent value  
 	* the solution is a sigular vector within V corresponding to the minimal  
 	* sigular value among S diagonal line coefficent， tipically x = V.col(11) as the repository.  
 	*  
 	* now we have got H, and we use QR decomposition to calculate K , R ,C  
 	*  
 	* H = [KR, -KRC] = [Ho, ho]  
 	* Ho = KR  
 	* ho = -KRC  
 	* so C = -Ho^-1 * ho  
 	* to our konwledge, in QR decomposition, R is a upper  
 	* triangular matrix, Q is a orthogonal matrix, and is our Ho =KR  
 	* matrix K is a upper triangular matrix, R is a orthogonal matrix  
 	*  
 	* so  we can make a samll trick H^-1 = R^-1 * K^-1 = QR  
 	* than we get rotation matrix R = Q^-1, and get K = R^-1  
 	* finnally we devide K by K(2,3) to get a normalized K.  
 	*  
**********/  
