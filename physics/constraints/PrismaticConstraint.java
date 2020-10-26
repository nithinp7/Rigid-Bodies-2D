
package physics.constraints;

import geometry.visual.Line;
import java.awt.Color;
import math.Vec2;
import physics.shapes.ConvPoly;
import static processing.core.PApplet.*;
import processing.core.PGraphics;

/**
 *
 * @author Nithin
 */
public final class PrismaticConstraint extends Constraint {

    private final Vec2 locAncA, locAncB, locNormA, locNormB;
    private final float dist;
    
    //globalDif must point from locAncA to locAncB
    public PrismaticConstraint(ConvPoly a, ConvPoly b, 
                               Vec2 locAncA, Vec2 locAncB,
                               Vec2 globalNormal) {
        
        super(a, b);
        
        this.locAncA = locAncA;
        this.locAncB = locAncB;
        
        dist = globalNormal.mag();
        globalNormal = globalNormal.normalize();
        
        locNormA = a.disorient(globalNormal);
        locNormB = b.disorient(globalNormal);
    }
    
    public PrismaticConstraint(ConvPoly a, ConvPoly b, Vec2 locAncA, Vec2 locAncB) {
        this(a, b, locAncA, locAncB, b.toGlobal(locAncB).sub(a.toGlobal(locAncA)));
    }
    
    @Override
    public void render(PGraphics g) {
        Line.drawLine(a.toGlobal(locAncA), b.toGlobal(locAncB), Color.GREEN, g);
    }
    
    @Override
    public float getBias(float dt) {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA));
        float errA = dif.dot(a.orient(locAncA)),
              errB = dif.dot(b.orient(locAncA));
        return -0.5f*abs(b.getAngle() - a.getAngle() - 0)/dt;
    }
    
    @Override
    public float[] getJacobian() {
        Vec2 dif = b.toGlobal(locAncB).sub(a.toGlobal(locAncA)),
             difN = dif.normalizeSafe(),
             normA = a.orient(locAncA),
             normB = b.orient(locAncB);
        float sepA = dif.dot(normA),
              sepB = -dif.dot(normB),
              perpA = normA.cross(dif),
              perpB = normB.cross(dif);
        return new float[] {
           locNormA.x, dif.y, dTheta,
           -dif.x, -dif.y, -dTheta
        };
    }
    
    @Override
    public float getLo() {
        return Float.NEGATIVE_INFINITY;
    }
    
    @Override
    public float getHi() {
        return Float.POSITIVE_INFINITY;
    }
}
