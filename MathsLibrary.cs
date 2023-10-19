using System.Collections;
using System.Collections.Generic;

public class MyVector3
{
    public float x, y, z;
    public MyVector3()
    {
        this.x = 0;
        this.y = 0;
        this.z = 0;
    }
    public MyVector3(float x, float y, float z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    public static Vector3 MyVector3ToVector3(MyVector3 vec)
    {
        Vector3 rv;
        rv.x = vec.x;
        rv.y = vec.y;
        rv.z = vec.z;
        return rv;
    }
    public static Vector3[] MyVector3ToVector3(MyVector3[] vec)
    {
        Vector3[] rv = new Vector3[vec.Length];
        for (int i = 0; i < vec.Length; i++)
        {
            rv[i].x = vec[i].x;
            rv[i].y = vec[i].y;
            rv[i].z = vec[i].z;
        }
        return rv;
    }
    public static MyVector3 Vector3ToMyVector3(Vector3 vec)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec.x;
        rv.y = vec.y;
        rv.z = vec.z;
        return rv;
    }
    public static MyVector3[] Vector3ToMyVector3(Vector3[] vec)
    {
        MyVector3[] rv = new MyVector3[vec.Length];
        for (int i = 0; i < vec.Length; i++)
        {
            rv[i] = new MyVector3(vec[i].x, vec[i].y, vec[i].z);
        }
        return rv;
    }
    public static MyVector3 AddVectors(MyVector3 vec1, MyVector3 vec2)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec1.x + vec2.x;
        rv.y = vec1.y + vec2.y;
        rv.z = vec1.z + vec2.z;
        return rv;
    }
    public static MyVector3 SubtractVectors(MyVector3 vec1, MyVector3 vec2)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec1.x - vec2.x;
        rv.y = vec1.y - vec2.y;
        rv.z = vec1.z - vec2.z;
        return rv;
    }
    public static MyVector3 ScaleVector(MyVector3 vec1, MyVector3 vec2)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec1.x * vec2.x;
        rv.y = vec1.y * vec2.y;
        rv.z = vec1.z * vec2.z;
        return rv;
    }
    public static MyVector3 DivideVector(MyVector3 vec1, float Scalar)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec1.x / Scalar;
        rv.y = vec1.y / Scalar;
        rv.z = vec1.z / Scalar;
        return rv;
    }
    public static MyVector3 operator +(MyVector3 ovec1, MyVector3 ovec2)
    {
        return AddVectors(ovec1, ovec2);
    }
    public static MyVector3 operator -(MyVector3 ovec1, MyVector3 ovec2)
    {
        return SubtractVectors(ovec1, ovec2);
    }
    public static MyVector3 operator *(MyVector3 ovec1, MyVector3 ovec2)
    {
        return ScaleVector(ovec1, ovec2);
    }
    public static MyVector3 operator *(MyVector3 ovec1, float scalar)
    {
        return ScaleVector(ovec1, new MyVector3(scalar, scalar, scalar));
    }
    public static MyVector3 operator *(float scalar, MyVector3 ovec1)
    {
        return ScaleVector(ovec1, new MyVector3(scalar, scalar, scalar));
    }
    public static MyVector3 operator /(MyVector3 ovec1, float divisor)
    {
        return DivideVector(ovec1, divisor);
    }
    public MyVector3 Normalize()
    {
        MyVector3 rv = new MyVector3();
        rv = this / this.GetLength();
        return rv;
    }
    public MyVector3 RotateVertexArountAxis(float angle, MyVector3 axis, MyVector3 vertex)
    {

        MyVector3 rv = (vertex * Mathf.Cos(angle)) + GetDotProduct(vertex, axis) * axis * (1 - Mathf.Cos(angle)) + CrossProduct(axis, vertex) * Mathf.Sin(angle);
        return rv;
    }
    public float GetLength()
    {
        float rv;
        rv = Mathf.Sqrt(x * x + y * y + z * z);
        return rv;
    }
    public float LengthSq()
    {
        float rv;
        rv = x * x + y * y + z * z;
        return rv;
    }

    public static float GetDotProduct(MyVector3 vec1, MyVector3 vec2, bool ShouldNormalize = true)
    {
        float rv;
        MyVector3 A = new MyVector3(vec1.x, vec1.y, vec1.z);
        MyVector3 B = new MyVector3(vec2.x, vec2.y, vec2.z);
        if (ShouldNormalize)
        {
            A = A.Normalize();
            B = B.Normalize();
        }
        rv = A.x * B.x + A.y * B.y + A.z * B.z;
        return rv;
    }
    public static MyVector3 EulerAngleToDirection(MyVector3 EulerAngles)
    {
        MyVector3 rv = new MyVector3();
        rv.x = Mathf.Cos(EulerAngles.y) * Mathf.Cos(EulerAngles.x);
        rv.y = Mathf.Sin(EulerAngles.x);
        rv.z = Mathf.Cos(EulerAngles.x) * Mathf.Sin(EulerAngles.y);
        return rv;
    }
    public static MyVector3 DirectionToEulerAngle(MyVector3 Direction)
    {
        MyVector3 rv = new MyVector3();
        rv.y = Mathf.Atan2(Direction.z, Direction.y);
        rv.y = Mathf.Atan2(Direction.x, Direction.z);
        rv.z = 0;
        return rv;
    }
    public static MyVector3 CrossProduct(MyVector3 vec1, MyVector3 vec2)
    {
        MyVector3 rv = new MyVector3();
        rv.x = vec1.y * vec2.z - vec1.z * vec2.y;
        rv.y = vec1.z * vec2.x - vec1.x * vec2.z;
        rv.z = vec1.x * vec2.y - vec1.y * vec2.x;
        return rv;
    }
    public static MyVector3 LerpVector(MyVector3 vec1, MyVector3 vec2, float t)
    {
        return vec1 * (1.0f - t) + vec2 * t;
    }
    public static MyVector3 FindCentre(MyVector3[] veclist)
    {
        MyVector3 rv = new MyVector3();
        float x = 0, y = 0, z = 0;
        for (int i = 0; i < veclist.Length; i++)
        {
            x += veclist[i].x;
            y += veclist[i].y;
            z += veclist[i].z;

            //Debug.Log(x + " " + y + " " + z);

        }
        rv = new MyVector3(x / veclist.Length, y / veclist.Length, z / veclist.Length);
        Debug.Log("av: " + rv.x + " " + rv.y + " " + rv.y + " length: " + veclist.Length);
        return rv;
    }
    
}
public class MyVector4
{
    public float x, y, z, w;
    public MyVector4()
    {
        this.x = 0;
        this.y = 0;
        this.z = 0;
        this.w = 0;
    }
    public MyVector4(float x, float y, float z, float w)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }
    
}
public class Matrix4by4 
{
    public float[,] values;
    public Matrix4by4()
    {
        values = new float[4, 4];
        this.values[0, 0] = 0;
        this.values[1, 0] = 0;
        this.values[2, 0] = 0;
        this.values[3, 0] = 0;

        this.values[0, 1] = 0;
        this.values[1, 1] = 0;
        this.values[2, 1] = 0;
        this.values[3, 1] = 0;

        this.values[0, 2] = 0;
        this.values[1, 2] = 0;
        this.values[2, 2] = 0;
        this.values[3, 2] = 0;

        this.values[0, 3] = 0;
        this.values[1, 3] = 0;
        this.values[2, 3] = 0;
        this.values[3, 3] = 0;
    }
    public static Matrix4by4 Identity
    {
        get
        {
            return new Matrix4by4(
                new MyVector4(1, 0, 0, 0),
                new MyVector4(0, 1, 0, 0),
                new MyVector4(0, 0, 1, 0),
                new MyVector4(0, 0, 0, 1)
                );
        }
    }
    public Matrix4by4(MyVector4 column1, MyVector4 column2, MyVector4 column3, MyVector4 column4)
    {
        values = new float[4, 4];
        values[0, 0] = column1.x;
        values[1, 0] = column1.y;
        values[2, 0] = column1.z;
        values[3, 0] = column1.w;

        values[0, 1] = column2.x;
        values[1, 1] = column2.y;
        values[2, 1] = column2.z;
        values[3, 1] = column2.w;

        values[0, 2] = column3.x;
        values[1, 2] = column3.y;
        values[2, 2] = column3.z;
        values[3, 2] = column3.w;

        values[0, 3] = column4.x;
        values[1, 3] = column4.y;
        values[2, 3] = column4.z;
        values[3, 3] = column4.w;
    }
    public Matrix4by4(MyVector3 column1, MyVector3 column2, MyVector3 column3, MyVector3 column4)
    {
        values = new float[4, 4];
        values[0, 0] = column1.x;
        values[1, 0] = column1.y;
        values[2, 0] = column1.z;
        values[3, 0] = 0;

        values[0, 1] = column2.x;
        values[1, 1] = column2.y;
        values[2, 1] = column2.z;
        values[3, 1] = 0;

        values[0, 2] = column3.x;
        values[1, 2] = column3.y;
        values[2, 2] = column3.z;
        values[3, 2] = 0;

        values[0, 3] = column4.x;
        values[1, 3] = column4.y;
        values[2, 3] = column4.z;
        values[3, 3] = 1;
    }
    public static MyVector4 operator *(Matrix4by4 matrix, MyVector4 scalar)
    {
        MyVector4 rv = new MyVector4();
        rv.x = matrix.values[0, 0] * scalar.x + matrix.values[0, 1] * scalar.y + matrix.values[0, 2] * scalar.z + matrix.values[0, 3] * scalar.w;
        rv.y = matrix.values[1, 0] * scalar.x + matrix.values[1, 1] * scalar.y + matrix.values[1, 2] * scalar.z + matrix.values[1, 3] * scalar.w;
        rv.z = matrix.values[2, 0] * scalar.x + matrix.values[2, 1] * scalar.y + matrix.values[2, 2] * scalar.z + matrix.values[2, 3] * scalar.w;
        rv.w = matrix.values[3, 0] * scalar.x + matrix.values[3, 1] * scalar.y + matrix.values[3, 2] * scalar.z + matrix.values[3, 3] * scalar.w;
        return rv;
    }
    public static MyVector3 operator *(Matrix4by4 matrix, MyVector3 scalar)
    {
        MyVector3 rv = new MyVector3();
        rv.x = matrix.values[0, 0] * scalar.x + matrix.values[0, 1] * scalar.y + matrix.values[0, 2] * scalar.z + matrix.values[0, 3] * 1;
        rv.y = matrix.values[1, 0] * scalar.x + matrix.values[1, 1] * scalar.y + matrix.values[1, 2] * scalar.z + matrix.values[1, 3] * 1;
        rv.z = matrix.values[2, 0] * scalar.x + matrix.values[2, 1] * scalar.y + matrix.values[2, 2] * scalar.z + matrix.values[2, 3] * 1;
        return rv;
    }
    public static Matrix4by4 operator *(Matrix4by4 m1, Matrix4by4 m2)
    {
        Matrix4by4 rm = new Matrix4by4();

        rm.values[0, 0] = m1.values[0, 0] * m2.values[0, 0] + m1.values[1, 0] * m2.values[0, 1] + m1.values[2, 0] * m2.values[0, 2] + m1.values[3, 0] * m2.values[0, 3];
        rm.values[1, 0] = m1.values[0, 0] * m2.values[1, 0] + m1.values[1, 0] * m2.values[1, 1] + m1.values[2, 0] * m2.values[1, 2] + m1.values[3, 0] * m2.values[1, 3];
        rm.values[2, 0] = m1.values[0, 0] * m2.values[2, 0] + m1.values[1, 0] * m2.values[2, 1] + m1.values[2, 0] * m2.values[2, 2] + m1.values[3, 0] * m2.values[2, 3];
        rm.values[3, 0] = m1.values[0, 0] * m2.values[3, 0] + m1.values[1, 0] * m2.values[3, 1] + m1.values[2, 0] * m2.values[3, 2] + m1.values[3, 0] * m2.values[3, 3];

        rm.values[0, 1] = m1.values[0, 1] * m2.values[0, 0] + m1.values[1, 1] * m2.values[0, 1] + m1.values[2, 1] * m2.values[0, 2] + m1.values[3, 1] * m2.values[0, 3];
        rm.values[1, 1] = m1.values[0, 1] * m2.values[1, 0] + m1.values[1, 1] * m2.values[1, 1] + m1.values[2, 1] * m2.values[1, 2] + m1.values[3, 1] * m2.values[1, 3];
        rm.values[2, 1] = m1.values[0, 1] * m2.values[2, 0] + m1.values[1, 1] * m2.values[2, 1] + m1.values[2, 1] * m2.values[2, 2] + m1.values[3, 1] * m2.values[2, 3];
        rm.values[3, 1] = m1.values[0, 1] * m2.values[3, 0] + m1.values[1, 1] * m2.values[3, 1] + m1.values[2, 1] * m2.values[3, 2] + m1.values[3, 1] * m2.values[3, 3];

        rm.values[0, 2] = m1.values[0, 2] * m2.values[0, 0] + m1.values[1, 2] * m2.values[0, 1] + m1.values[2, 2] * m2.values[0, 2] + m1.values[3, 2] * m2.values[0, 3];
        rm.values[1, 2] = m1.values[0, 2] * m2.values[1, 0] + m1.values[1, 2] * m2.values[1, 1] + m1.values[2, 2] * m2.values[1, 2] + m1.values[3, 2] * m2.values[1, 3];
        rm.values[2, 2] = m1.values[0, 2] * m2.values[2, 0] + m1.values[1, 2] * m2.values[2, 1] + m1.values[2, 2] * m2.values[2, 2] + m1.values[3, 2] * m2.values[2, 3];
        rm.values[3, 2] = m1.values[0, 2] * m2.values[3, 0] + m1.values[1, 2] * m2.values[3, 1] + m1.values[2, 2] * m2.values[3, 2] + m1.values[3, 2] * m2.values[3, 3];

        rm.values[0, 3] = m1.values[0, 3] * m2.values[0, 0] + m1.values[1, 3] * m2.values[0, 1] + m1.values[2, 3] * m2.values[0, 2] + m1.values[3, 3] * m2.values[0, 3];
        rm.values[1, 3] = m1.values[0, 3] * m2.values[1, 0] + m1.values[1, 3] * m2.values[1, 1] + m1.values[2, 3] * m2.values[1, 2] + m1.values[3, 3] * m2.values[1, 3];
        rm.values[2, 3] = m1.values[0, 3] * m2.values[2, 0] + m1.values[1, 3] * m2.values[2, 1] + m1.values[2, 3] * m2.values[2, 2] + m1.values[3, 3] * m2.values[2, 3];
        rm.values[3, 3] = m1.values[0, 3] * m2.values[3, 0] + m1.values[1, 3] * m2.values[3, 1] + m1.values[2, 3] * m2.values[3, 2] + m1.values[3, 3] * m2.values[3, 3];
        return rm;
    }
    public static Matrix4by4 ScaleMatrix(float x, float y, float z)
    {
        Matrix4by4 scalar;
        scalar = new Matrix4by4(new MyVector3(x,0,0), new MyVector3(0,y,0), new MyVector3(0,0,z), new MyVector3());
        return scalar;
    }
    public static Matrix4by4 TranslateMatrix(float x, float y, float z)
    {
        Matrix4by4 translator;
        translator = new Matrix4by4(new MyVector3(1, 0, 0), new MyVector3(0, 1, 0), new MyVector3(0, 0, 1), new MyVector3(x, y, z));
        return translator;
    }
    public static Matrix4by4 YawMatrix(float angle)
    {
        Matrix4by4 yawmatrix = new Matrix4by4(
            new MyVector3(Mathf.Cos(angle), 0, -Mathf.Sin(angle)),
            new MyVector3(0, 1, 0),
            new MyVector3(Mathf.Sin(angle), 0, Mathf.Cos(angle)),
            new MyVector3()
            );
        return yawmatrix;
    }
    public static Matrix4by4 PitchMatrix(float angle)
    {
        Matrix4by4 pitchmatrix = new Matrix4by4(
            new MyVector3(1,0,0),
            new MyVector3(0,Mathf.Cos(angle),Mathf.Sin(angle)),
            new MyVector3(0,-Mathf.Sin(angle), Mathf.Cos(angle)),
            new MyVector3()
            );
        return pitchmatrix;
    }
    public static Matrix4by4 RollMatrix(float angle)
    {
        Matrix4by4 rollmatrix = new Matrix4by4(
            new MyVector3(Mathf.Cos(angle), Mathf.Sin(angle),0),
            new MyVector3(-Mathf.Sin(angle), Mathf.Cos(angle), 0),
            new MyVector3(0,0,1),
            new MyVector3()
            );
        return rollmatrix;
    }
    public static Matrix4by4 InverseRotationMatrix(Matrix4by4 rotator)
    {

        Matrix4by4 rv = new Matrix4by4();
        rv.values[0, 0] = rotator.values[0, 0];
        rv.values[0, 1] = rotator.values[1, 0];
        rv.values[0, 2] = rotator.values[2, 0];
        rv.values[0, 3] = rotator.values[3, 0];

        rv.values[1, 0] = rotator.values[0, 1];
        rv.values[1, 1] = rotator.values[1, 1];
        rv.values[1, 2] = rotator.values[2, 1];
        rv.values[1, 3] = rotator.values[3, 1];

        rv.values[2, 0] = rotator.values[0, 2];
        rv.values[2, 1] = rotator.values[1, 2];
        rv.values[2, 2] = rotator.values[2, 2];
        rv.values[2, 3] = rotator.values[3, 2];

        rv.values[3, 0] = rotator.values[0, 3];
        rv.values[3, 1] = rotator.values[1, 3];
        rv.values[3, 2] = rotator.values[2, 3];
        rv.values[3, 3] = rotator.values[3, 3];
        return rv;
    }
    public static Matrix4by4 InverseTranslationMatrix(Matrix4by4 translator)
    {
        Matrix4by4 rv = translator;
        rv.values[0, 3] *= -1;
        rv.values[1, 3] *= -1;
        rv.values[2, 3] *= -1;
        return rv;
    }
    public static Matrix4by4 InverseScaleMatrix(Matrix4by4 scalar)
    {
        Matrix4by4 rv = scalar;
        rv.values[0,0] = 1 / rv.values[0,0];
        rv.values[1,1] = 1 / rv.values[1,1];
        rv.values[2,2] = 1/ rv.values[2,2];
        return rv;
    }
}
public class Quat
{
    public float x, y, z, w;
    public Quat()
    {
        x = 0;
        x = 0;
        y = 0;
        w = 0;
    }
    public Quat(float angle, MyVector3 axis)
    {
        float halfangle = angle / 2;
        w = Mathf.Cos(halfangle);
        x = axis.x * Mathf.Sin(halfangle);
        y = axis.y * Mathf.Sin(halfangle);
        z = axis.z * Mathf.Sin(halfangle);
    }
    public static Quat operator *(Quat q1, Quat q2)
    {
        Quat rq = new Quat();
        rq.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
        rq.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z + q1.z * q2.y;
        rq.y = q1.w * q2.y + q1.y * q2.w + q1.z + q2.x + q1.x * q2.z;
        rq.z = q1.w * q2.z + q1.z * q2.w + q1.x + q2.y + q1.y * q2.x;
        return rq;
    }
    public MyVector4 GetAxisAngle()
    {
        MyVector4 rv = new MyVector4();
        float halfangle = Mathf.Acos(w);
        rv.w = halfangle * 2;

        rv.x = y / Mathf.Sin(halfangle);
        rv.y = y / Mathf.Sin(halfangle);
        rv.z = z / Mathf.Sin(halfangle);

        return rv;
    }
    public Quat Inverse()
    {
        Quat rq = new Quat();
        float num = (x * x) + (y * y) + (z * z) + (w * w);
        rq.x = -x * 1 / num;
        rq.y = -y * 1 / num;
        rq.z = -z * 1 / num;
        rq.w = w * 1 / num;
        return rq;
    }
    public static Quat SLERP(Quat q, Quat r, float t)
    {
        t = Mathf.Clamp(t, 0.01f, 1f);
        Quat d = r * q.Inverse();
        MyVector4 axisangle = d.GetAxisAngle();
        Quat dT = new Quat(axisangle.w * t, new MyVector3(axisangle.x, axisangle.y, axisangle.z));
        return dT * q;
    }
    public static Quat Euler(float z, float y, float x)
    {
        Quat rq = new Quat();
        rq.x = Mathf.Sin(z / 2) * Mathf.Cos(y / 2) * Mathf.Cos(x / 2) - Mathf.Cos(z / 2) * Mathf.Sin(y / 2) * Mathf.Sin(x / 2);
        rq.y = Mathf.Cos(z / 2) * Mathf.Sin(y / 2) * Mathf.Cos(x / 2) + Mathf.Sin(z / 2) * Mathf.Cos(y / 2) * Mathf.Sin(x / 2);
        rq.z = Mathf.Cos(z / 2) * Mathf.Cos(y / 2) * Mathf.Sin(x / 2) - Mathf.Sin(z / 2) * Mathf.Sin(y / 2) * Mathf.Cos(x / 2);
        rq.w = Mathf.Cos(z / 2) * Mathf.Cos(y / 2) * Mathf.Cos(x / 2) + Mathf.Sin(z / 2) * Mathf.Sin(y / 2) * Mathf.Sin(x / 2);

        return rq;
    }
    public static Quaternion QuatToQuaternion(Quat q)
    {
        Quaternion rq = new Quaternion();
        rq.w = q.w;
        rq.x = q.x;
        rq.y = q.y;
        rq.z = q.z;
        return rq;
    }
}
public class AABB
{
    MyVector3 MinExtent;
    MyVector3 MaxExtent;
    public AABB(MyVector3 min, MyVector3 max)
    {
        MinExtent = min;
        MaxExtent = max;

    }
    public float Top
    {
        get { return MaxExtent.y; }
    }
    public float Bottom
    {
        get { return MinExtent.y; }
    }
    public float Left
    {
        get { return MinExtent.x; }
    }
    public float Right
    {
        get { return MaxExtent.x; }
    }
    public float Front
    {
        get { return MaxExtent.z; }
    }
    public float Back
    {
        get { return MinExtent.z; }
    }
    public static bool BoxIntersection(AABB box1, AABB box2)
    {
        return !(box2.Left > box1.Right
            || box2.Right < box1.Left
            || box2.Top < box1.Bottom
            || box2.Bottom > box1.Top
            || box2.Back > box1.Front
            || box2.Front < box1.Back);
    }
    public static bool LineIntersection(AABB box, MyVector3 startpoint, MyVector3 endpoint, out MyVector3 intersectionpoint)
    {
        float lowest = 0.0f;
        float highest = 1.0f;
        intersectionpoint = new MyVector3();
        if (!IntersectingAxis('x', box, startpoint, endpoint, ref lowest, ref highest)) return false;
        if (!IntersectingAxis('y', box, startpoint, endpoint, ref lowest, ref highest)) return false;
        if (!IntersectingAxis('z', box, startpoint, endpoint, ref lowest, ref highest)) return false;
        intersectionpoint = MyVector3.LerpVector(startpoint, endpoint, lowest);
        return true;
    }
    public static bool IntersectingAxis(char axis, AABB box, MyVector3 startpoint, MyVector3 endpoint, ref float lowest, ref float highest)
    {
        float minimum = 0.0f, maximum = 0.0f;
        if (axis == 'x')
        {
            minimum = (box.Left - startpoint.x) / (endpoint.x - startpoint.x);
            maximum = (box.Right - startpoint.x) / (endpoint.x - startpoint.x);
        }
        else if (axis == 'y')
        {
            minimum = (box.Bottom - startpoint.y) / (endpoint.y - startpoint.y);
            maximum = (box.Top - startpoint.y) / (endpoint.y - startpoint.y);
        }
        else if (axis == 'z')
        {
            minimum = (box.Back - startpoint.z) / (endpoint.z - startpoint.z);
            maximum = (box.Front - startpoint.z) / (endpoint.z - startpoint.z);
        }
        if (maximum < minimum)
        {
            float temp = maximum;
            maximum = minimum;
            minimum = temp;
        }

        if (maximum < lowest) return false;
        if (minimum > highest) return false;

        lowest = Mathf.Max(minimum, lowest);
        highest = Mathf.Min(maximum, highest);

        if (lowest > highest) return false;

        return true;
    }
}
public class BoundingSphere
{
    public MyVector3 centrepoint;
    public float radius;
    public BoundingSphere(MyVector3 Centre, float Radius)
    {
        this.centrepoint = Centre;
        this.radius = Radius;
    }
    public static float SqDistanceFromPointToSegment(MyVector3 a, MyVector3 b, MyVector3 c)
    {
        float distancesq = 0;
        MyVector3 ab = b - a;
        MyVector3 ac = c - a;
        distancesq  = ac.LengthSq() - (MyVector3.GetDotProduct(ac,ab,false) * MyVector3.GetDotProduct(ac, ab,false)/ab.LengthSq());
        return distancesq;
    }
    public static bool LineIntersection(BoundingSphere sphere, MyVector3 start, MyVector3 end, out MyVector3 intersect, out float z)
    {
        float distancesq = SqDistanceFromPointToSegment(start, end, sphere.centrepoint);
        float intersectpoint = MyVector3.GetDotProduct(sphere.centrepoint-start, end-start, false) / (end-start).LengthSq();

        //equation for point on circle y = sqr r(r^2 - x^2)
        //z will always be the same at a given x
        z = Mathf.Sqrt((sphere.radius * sphere.radius) - distancesq);
        intersect = (start + end) * intersectpoint;
        return distancesq < sphere.radius * sphere.radius;
    }
}
