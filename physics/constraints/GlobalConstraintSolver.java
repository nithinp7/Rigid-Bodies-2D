
package physics.constraints;

import math.Mat;
import math.Vec;

/**
 *
 * @author Nithin
 */
public class GlobalConstraintSolver {
    
    private static final int iters=20;
    
    private Mat M_INV=new Mat(0,0);
    
    public GlobalConstraintSolver() {}
    
    public void setMassInvMatrix(float[] massInvEntries) {
        M_INV = Mat.diagonal(massInvEntries);
    }
    
    //returns constraint impulses
    //TODO:
    //generalize to arbitrary hi and lo bounds on x 
    public Vec solveConstraintImpulses(Mat J, Vec V, Vec dV, Vec bias,
                                       int[] limDep, float[] lo, float[] hi, 
                                       float cfm) {
        Vec V_t = V.add(dV);
        
        Mat Jt = J.transpose();
        Mat A = J.mul(M_INV).mul(Jt);
        Vec b = J.mul(V_t).add(bias).mul(-1);
        
        float[] x = new float[A.rows];
        
        //NOW WE SOLVE Ax = b
        
        for(int k=0; k<iters; k++) 
            for(int i=0; i<A.rows; i++) {
                float delta = 0;
                for(int j=0; j<i; j++)
                    delta += A.get(i,j)*x[j];
                for(int j=i+1; j<A.rows; j++) 
                    delta += A.get(i,j)*x[j];
                x[i] = (b.get(i)-delta)/(A.get(i,i)+cfm);
                float s=1;
                if(limDep[i]>=0)
                    s=x[limDep[i]];
                if(x[i]<lo[i]*s)
                    x[i]=lo[i]*s;
                if(x[i]>hi[i]*s)
                    x[i]=hi[i]*s;
            }
        
        Vec Pc = Jt.mul(new Vec(x));
        
        return Pc;//V_t.add(M_INV.mul(Pc));
    }
}
