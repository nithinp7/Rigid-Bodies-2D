
package physics.constraints;

import geometry.collision.Manifold;
import physics.shapes.ConvPoly;
import java.util.ArrayList;
import processing.core.PGraphics;
import static geometry.collision.CollisionDetector.*;
import math.*;
import static processing.core.PApplet.abs;

/**
 *
 * @author Nithin
 */

public final class ConstraintsArbiter {
    
    /**
     * NOTE:
     * non-penetration and friction constraint jacobians are retrieved from
     * Manifold objects, all other types of constraint jacobians are retrieved
     * from Constraint objects
     */
    
    private final ArrayList<ConvPoly> polys = new ArrayList<>();
    private final ArrayList<Manifold> manifolds = new ArrayList<>();
    private final ArrayList<Constraint> constraints = new ArrayList<>();
    
    private final GlobalConstraintSolver solver;
    
    //TODO: decide if CollisionDetector should only be used statically
    //private final CollisionDetector detector;
    
    public class ColliderIndex {
        
        private final int index; 
        
        private ColliderIndex(int i) {
            index = i;
        }
    }
    
    public ConstraintsArbiter() {
        solver = new GlobalConstraintSolver();
    }
    
    public void removeManifold(ConvPoly a, ConvPoly b) {
        manifolds.removeIf(m -> m.getA()==a && m.getB()==b ||
                                m.getA()==b && m.getB()==a);
    }
    
    public void addPoly(ConvPoly newPoly) {
        newPoly.colIndex = new ColliderIndex(polys.size());
        polys.stream().forEach(p -> manifolds.add(new Manifold(p,newPoly)));
        polys.add(newPoly);
        
        float[] M_INV = new float[3*polys.size()];
        for(int i=0; i<polys.size(); i++) {
            ConvPoly p = polys.get(i);
            M_INV[3*i] = p.massInv;
            M_INV[3*i+1] = p.massInv;
            M_INV[3*i+2] = p.momInertiaInv;
        }
        solver.setMassInvMatrix(M_INV);
    }
    
    public void addConstraint(Constraint c) {
        constraints.add(c);
    }
    
    public void clear() {
        polys.clear();
        manifolds.clear();
        constraints.clear();
    }
    
    public void update(float dt) {
        manifolds.stream().forEach(m -> checkSAT(m));
        
        Manifold[] C = manifolds
                         .stream()
                         .filter(Manifold::isColliding)
                         .toArray(Manifold[]::new);
        
        int n = polys.size(),
            cLen = manifolds
                    .stream()
                    .map(Manifold::getContactCount)
                    .reduce((Integer a, Integer b) -> a+b)
                    .orElse(0),
            contactsCount = 2*cLen + constraints.size();
        
        if(contactsCount==0) return;
        
        float[] bias = new float[contactsCount],
                lo =  new float[contactsCount],
                hi =  new float[contactsCount];
        int[] limDep = new int[contactsCount];
        
        float[][] J = new float[contactsCount][n*3];
        
        //CONTACT AND FRICTION CONSTRAINTS
        for(int i=0,offs=0; i<cLen; i++) {
            Manifold m = C[i-offs];
            ConvPoly a=m.getRef(), b=m.getInc();
            int ai=a.colIndex.index, bi=b.colIndex.index;
            float[][] js = m.getContactsJacobian();
            
            //iterate through all contact points on this manifold
            for(int cp=0; cp<js.length; cp++) {
                i += cp;
                offs += cp;
                
                float[] j = js[cp];
                
                J[i][3*ai] = j[0];
                J[i][3*ai+1] = j[1];
                J[i][3*ai+2] = j[2];

                J[i][3*bi] = j[3];
                J[i][3*bi+1] = j[4];
                J[i][3*bi+2] = j[5];

                //slop
                float d = abs(j[12])<0.1f?0:j[12];
                bias[i] = m.getRestitution()*m.getNormalVelocity() - 0.2f*d/dt;//m.getRestitution()*m.getNormalVelocity()
                lo[i] = 0;
                hi[i] = Float.POSITIVE_INFINITY;
                limDep[i] = -1;

                J[i+cLen][3*ai] = j[6];
                J[i+cLen][3*ai+1] = j[7];
                J[i+cLen][3*ai+2] = j[8];

                J[i+cLen][3*bi] = j[9];
                J[i+cLen][3*bi+1] = j[10];
                J[i+cLen][3*bi+2] = j[11];

                bias[i+cLen] = 0;
                lo[i+cLen] = -0.6f/js.length;
                hi[i+cLen] = 0.6f/js.length;
                limDep[i+cLen] = i;
            }
        }
        
        for(int i=0; i<constraints.size(); i++) {
            Constraint c = constraints.get(i);
            int ai=c.a.colIndex.index, bi=c.b.colIndex.index;
            float[] j = c.getJacobian();
            
            if(j==null) continue;
            
            J[i+2*cLen][3*ai] = j[0];
            J[i+2*cLen][3*ai+1] = j[1];
            J[i+2*cLen][3*ai+2] = j[2];
            
            J[i+2*cLen][3*bi] = j[3];
            J[i+2*cLen][3*bi+1] = j[4];
            J[i+2*cLen][3*bi+2] = j[5];
            
            bias[i+2*cLen] = c.getBias(dt);
            lo[i+2*cLen] = c.getLo();
            hi[i+2*cLen] = c.getHi();
            limDep[i+2*cLen] = -1;
        }
        
        Mat Jmat = new Mat(contactsCount, n*3,J);
        
        float[] V = new float[n*3],
                dV = new float[n*3];
        
        for(int i=0; i<n; i++) {
            ConvPoly p = polys.get(i);
            Vec2 v = p.getVelocity(), 
                 dv = p.getDeltaVelocity(),
                 f = p.getForceExt();
            float omega = p.getAngularVelocity(),
                  domega = p.getDeltaAngularVelocity(),
                  torque = p.getTorqueExt();
            
            V[3*i] = v.x;
            V[3*i+1] = v.y;
            V[3*i+2] = omega;
            
            dV[3*i] = f.x*p.massInv*dt + dv.x;
            dV[3*i+1] = f.y*p.massInv*dt + dv.y;
            dV[3*i+2] = torque*p.momInertiaInv + domega;
        }
        
        Vec Pc = solver.solveConstraintImpulses(
                                Jmat, 
                                new Vec(V), 
                                new Vec(dV),
                                new Vec(bias),
                                limDep,
                                lo,
                                hi,
                                0.0001f
                );
        
        for(int i=0; i<n; i++) {
            ConvPoly p = polys.get(i);
            p.applyImpulse(Pc.get(3*i),
                           Pc.get(3*i+1),
                           Pc.get(3*i+2));
        }
    }
    
    //TODO: filter manifolds more thoroughly before deciding to render;
    //allow for dormant manifolds
    //NOTE: for debug visualisations only
    public void render(PGraphics g) {
        manifolds.stream()
                 .filter(m -> m.isActive() && m.isColliding())
                 .forEach(m -> m.render(g));
        constraints.stream()
                   .forEach(c -> c.render(g));
    }
    
}
