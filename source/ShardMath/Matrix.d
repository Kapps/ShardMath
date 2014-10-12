module ShardMath.Matrix;
private import std.path;
private import std.traits;
import std.conv;
import ShardMath.Vector;
import std.exception;
import std.math;
import std.c.string;
import std.format;
import tested;

@name("ToString Tests")
unittest {			
	Matrix2f test = Matrix2f(1, 2, 3, 4);
	string expected = "[1, 2]\n[3, 4]";
	assert(test.text == expected);
}

@name("Constructor Tests")
unittest {
	Matrix2f test1 = Matrix2f(1, 2, 3, 4);
	assert(test1.m11 == 1);
	assert(test1.m12 == 2);
	assert(test1.m21 == 3);
	assert(test1.m22 == 4);
	auto t2 = Matrix3f(cast(float[3][3])[cast(float[3])[1f, 2f, 3f], cast(float[3])[4f, 5f, 6f], cast(float[3])[7f, 8f, 9f]]);
	assert(t2.m11 == 1);
	assert(t2.m21 == 4);
	assert(t2.m33 == 9);
	auto t3 = Matrix2f(1, 1, 1, 1);
	assert(t3.m11 == 1);
	assert(t3.m12 == 1);
	assert(t3.m21 == 1);
}

@name("Trace Tests")
unittest {
	auto mat = Matrix2f(1, 2, 3, 4);
	assert(mat.trace == 5);
	assert(Matrix4f.identity.trace == 4);
}

@name("Determinant Tests")
unittest {
	assert(Matrix2f(
		-5, 1, 
		1, 3
		).determinant == -16);
	float det = Matrix2f(
		4, -1, 
		2, -1/2f
		).determinant;
	assert(abs(det) < 0.0001f);
}

@nogc:

/// Represents a row-major square NxN coefficient matrix with helper methods for graphics programming.
struct Matrix(int N, T) if(N >= 2) {

	/// Implements opApply, either by iterating over all values alone, or by iterating over (row, column, value).
	int opApply(int delegate(ref T) dg) {
		int result;
		foreach(ref T element; this.elementsSingleDim)
			if((result = dg(element)) != 0)
				break;
		return result;
	}
	
	/// Ditto
	int opApply(int delegate(size_t, size_t, ref T) dg) {
		int result;
		size_t index;
		for(size_t row = 0; row < N; row++) {
			for(size_t col = 0; col < N; col++) {
				if((result = dg(row, col, elementsSingleDim[index++])) != 0)
					break;
			}
		}
		return result;
	}

	/// Returns this vector as a multi-line string in the same way that an array would be.
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const {
		for(int y = 0; y < N; y++) {
			sink("[");
			for(int x = 0; x < N; x++) {
				sink.formatValue(elementsSingleDim[y * N + x], fmt);
				if(x != N - 1)
					sink(", ");
			}
			sink("]");
			if(y != N - 1)
				sink("\n");
		}
	}

@safe pure nothrow:

	/// If T is a floating point type, this is alased to T; otherwise, float.
	/// This is used for operations that do not have a meaningful result for integer matrices.
	alias Select!(isFloatingPoint!(T), T, float) DecimalType;	

	// For easier external access.
	alias N NumRows;
	alias N NumColumns;
	alias T ElementType;

	/// Creates a matrix from an already created static array.
	this(in T[N][N] elements) { 
		this.elements = elements; 
	}

	/// Creates a matrix where all of the major elements have a given value.
	this(T value) { 
		for(int i = 0; i < N; i++) 
			this.elementsSingleDim[i] = value; 
	}

	/// Creates a matrix by specifying each value directly (m11 -> mNN).
	this(T...)(T elements) if(T.length == N * N) {
		foreach(i, val; elements)
			elementsSingleDim[i] = val;
	}

	/// Gets an instance of the identity matrix.
	@property static Matrix!(N, T) identity() @trusted {
		auto res = _identity;
		return res;
	}

	// TODO: Try to implement pass-by-ref.

	/// Binary operators for scalars and matrices.
	Matrix opBinary(string op)(in Matrix second) if(op == "+" || op == "-" || op == "*" || op == "/") {
		Matrix result;
		mixin(GenOp(op, "result", "this", "second"));
		return result;
	}

	/// Ditto
	void opBinaryAssign(string op)(in Matrix other) if(op == "+" || op == "-" || op == "*" || op == "/") {		
		Matrix result = this.opBinary!op(other);
		this = result;
		return this;
		// TODO: Bring this back!
		//mixin(GenOpAssign(op, "this", "other"));
		//return this;
	}

	/// Implementations for basic operators.
	T opIndex(size_t row, size_t column) const {
		return this.elementsSingleDim[row * N + column];
	}

	/// Ditto
	void opIndexAssign(size_t row, size_t column, T value) {
		this.elementsSingleDim[row * N + column] = value;
	}

	/// Ditto
	bool opEquals(in Matrix other) {
		for(size_t i = 0; i < elementsSingleDim.length; i++)
			if(!approxEqual(elementsSingleDim[i], other.elementsSingleDim[i]))
				return false;
		return true;
	}

	/// Ditto
	void opAssign(in Matrix other) @trusted {
		memcpy(this.elementsSingleDim.ptr, other.elementsSingleDim.ptr, T.sizeof * N * N);
	}

	static if(N == 4) {		
		/// Gets a Vector pointing in the given direction after rotations are applied by this Matrix.
		@property Vector!(3, T) right() const {
			return Vector!(3, T)(m11, m12, m13);	
		}
		
		/// Ditto
		@property Vector!(3, T) left() const {
			return -right;
		}

		/// Ditto
		@property Vector!(3, T) forward() const {
			return -backward;
		}

		/// Ditto
		@property Vector!(3, T) backward() const {
			return Vector!(3, T)(m31, m32, m33);
		}

		/// Gets or sets the position component of this matrix.		
		@property Vector!(3, T) translation() const {
			return Vector!(3, T)(m41, m42, m43);
		}

		/// Ditto
		@property void translation(in Vector!(3, T) value) {
			m41 = value.x;
			m42 = value.y;
			m43 = value.z;		
		}
	}

	/// Returns a transposed version of this matrix.
	@property Matrix!(N, T) transposed() const {
		Matrix result = void;		
		mixin(transposeMixin());
		return result;
	}

	private static string transposeMixin() {
		string MixinString = "";
		for(int y = 1; y <= N; y++) {
			for(int x = 1; x <= N; x++) {				
				MixinString ~= ("result.m" ~ to!string(y) ~ to!string(x) ~ " = this.m" ~ to!string(x) ~ to!string(y) ~ ";");
			}
		}	
		return MixinString;	
	}

	/+
	/// Calculates the Minor of this Matrix for the given row and column.
	/// That is, this Matrix with the given row and column removed.
	/// Params:
	/// 	row = The row to calculate the cofactor for.
	/// 	column = The column to calculate the cofactor for.
	T Minor(int row, int column) {
		assert(row <= Numrows && column <= Numcolumns && row > 0 && column > 0);
		alias Matrix!(N - 1, T) RetType;
		RetType result = void;
		for(int y = 1; y <= N; y++) {
			for(int x = 1; x <= N; x++) {
				
			}
		}
		mixin(MinorMixin(row, column));
		return result.determinant;		
	}

	private static string MinorMixin(int row, int column) {
		static assert(0, "Change to using elements instead of mixing in stuff we don't even know at compile-time.");
		// Remember, Minor = det(Matrix without row, column), Cofactor = -1^(row+column) * Minor(row, column).
		string result = "int targetrow = 0, targetcol = 0; ";
		for(int Y = 1; Y < N; Y++) {
			for(int X = 1; X < N; X++) {
				result ~= "targetrow = Y >= row ? Y - 1 : Y; ";
				result ~= "targetcol = X >= column ? X - 1 : X; ";
				result ~= "result.m" ~ to!string(targetrow) ~ to!string(targetcol) ~ " = this.m" ~ to!string(Y) ~ to!string(X) ~ "; ";
			}
		}
		return result;
	}	

	/// Calculates the Cofactor of this Matrix for the given row and column.
	/// Params:
	/// 	row = The row to calculate the cofactor for.
	/// 	column = The column to calculate the cofactor for.
	T Cofactor(int row, int column) {
		return (row + column) % 2 == 0 ? Minor(row, column) : -Minor(row, column);
	}

	unittest {
		Matrix4f testMat = Matrix4f(
			 4,  0, 10,  4,
			-1,  2,  3,  9,
			 5, -5, -1,  6,
			 3,	 7,  1, -2
		);
		Matrix3f Actual = testMat.minor(3, 3);
		Matrix3f expected = Matrix3f(
			 4,  0,  4,
			-1,  2,  9,
			 3,  7, -2
		);
		assert(Actual == expected);
		assert(testMat.Cofactor(3, 3) == Actual);
	}+/
	
	/// Calculates the determinant of this matrix.
	@property T determinant() {		
		static if(N == 2) {
			return (m11 * m22) - (m12 * m21);			
		} else static if(N == 3) {
			return (m11 * m22 * m33) + (m12 * m23 * m31) + (m13 * m21 * m32) - (m12 * m21 * m33) - (m11 * m23 * m32) - (m13 * m22 * m31);
		} else static if(N == 4) {
			assert(0, "Not yet implemented.");
			// Calculate with method of cofactors.
		} else {
			assert(0, "Not yet implemented.");
		}		
	}

	static if(isFloatingPoint!(T)) {
		// TODO: Support DecimalType. Keep in mind that InvertInline wouldn't work.

		version(None) { // TODO: Broked.
		/// Calculates the inverse of this Matrix.
		@property Matrix!(N, T) Inverse() const {
			// TODO: Consider optimizing.
			Matrix!(N, T) result = this;			
			result.InvertInline();
			return result;			
		}

		unittest {
			Matrix2f Inverted = Matrix2f(-5, 1, 1, 3).InvertInline();
			assert(Inverted == Matrix2f(-3/16f, 1/16f, 1/16f, 5/16f));
		}

		/// Inverts this Matrix by altering it's own elements, returning the same instance of the Matrix.
		ref Matrix InvertInline() {
			assert(IsInvertible);
			// TODO: Can optimize this by storing the results of the calculations we used for determinant.
			static if(N == 2) {				
				T Det = determinant;
				assert(Det != 0);
				T InvDet = 1f / Det;			
				m22 *= InvDet;
				m12 *= -InvDet;
				m21 *= -InvDet;
				m11 *= InvDet;
				return this;
			}
			assert(0, "Not yet implemented.");
		}
		}
	}

	/// Indiciates whether this Matrix can be inverted.
	@property bool isInvertible() {
		return abs(determinant) > 0.0001f;
	}
	 
	private static string elementMixin() {
		string result = "union { T[N * N] elementsSingleDim; T[N][N] elements;\r\nstruct {";
		for(int y = 1; y <= N; y++) {
			result ~= "union {\r\n Vector!(N, T) row" ~ to!string(y) ~ ";\r\nstruct {";
			for(int x = 1; x <= N; x++) {
				string element = "m" ~ to!string(y) ~ to!string(x);
				result ~= "T " ~ element ~ ";\r\n";
			}
			result ~= "}\r\n}";
		}
		result ~= "}\r\n}";
		return result;
	}

	mixin(elementMixin());	

private:
	__gshared immutable Matrix!(N, T) _identity;	
	
	shared static this() {		
		for(size_t row = 0; row < N; row++) {
			for(size_t col = 0; col < N; col++) {
				if(row == col)
					_identity.elements[row][col] = 1;
				else
					_identity.elements[row][col] = 0;
			}
		}
		for(size_t i = 0; i < N; i++) {
			_identity.elements[i][i] = 1;
		}
	}	

	static string GenOp(string operator, string appliedTo, string left, string right) {
		// TODO: SIMD		
		//static assert(operator.length == 1, "expected single operator, such as + or -.");
		string result = "";		
		for(int x = 1; x <= N; x++) {
			for(int y = 1; y <= N; y++) {
				string element = "m" ~ to!string(y) ~ to!string(x);
				if(operator != "*")
					result ~= appliedTo ~ "." ~ element ~ " = " ~ left ~ "." ~ element ~ " " ~ operator ~ " " ~ right ~ "." ~ element ~ "; ";
				else {
					result ~= appliedTo ~ "." ~ element ~ " = 0; ";
					for(int z = 1; z <= N; z++) {
						string leftElement = "m" ~ to!string(x) ~ to!string(z);
						string rightElement = "m" ~ to!string(z) ~ to!string(y);
						result ~= appliedTo ~ "." ~ element ~ "+= " ~ left ~ "." ~ leftElement ~ " * " ~ right ~ "." ~ rightElement ~ "; ";
					}
				}					
			}
		}		
		return result;
	}	

	static string GenOpAssign(string operator, string left, string right) {
		string leftAccess = left == "this" ? "" : left ~ ".";
		// TODO: SIMD		
		//static assert(operator.length == 1, "expected single operator, such as + or -.");
		string result = "";		
		// TODO: Less hackish resulting in actual performance benefits for inlining, not costs...
		if(operator == "*") {			
			enum string tmpMatrixName = "__tmp_MatrixNT_opAssign";
			if(leftAccess != "")
				result ~= "Matrix " ~ tmpMatrixName ~ " = *" ~ leftAccess ~ " * " ~ right ~ "; ";
			else
				result ~= "Matrix " ~ tmpMatrixName ~ " = Matrix.opBinary!(\"*\")(" ~ right ~ "); ";
		} else {
			for(int x = 1; x <= N; x++) {
				for(int y = 1; y <= N; y++) {								
					string element = "m" ~ to!string(x) ~ to!string(y);					
					if(operator != "*")
						result ~= leftAccess ~ element ~ " " ~ operator ~ "= " ~ right ~ "." ~ element ~ "; ";
					else
						result ~= leftAccess ~ element ~ " = " ~ right ~ "." ~ element ~ "; ";
				}
			}		
		}				
		return result;
	}		
	
	/// Calculates the trace of this Matrix, usually noted by tr(Matrix).
	@property T trace() {
		T sum = 0;
		for(size_t i = 0; i < N * N; i += N + 1)
			sum += elementsSingleDim[i];
		return sum;
	}
}

/// creates a perspective field of view matrix with the given parameters.
/// This Matrix is valid as both a DirectX and OpenGL Projection Matrix, but for OpenGL must be transposed.
/// Note that glUniformMatrix4fv contains a transpose parameter.
/// Also note that using a transposed version of this for OpenGL results in different GLSL multiplication order (Model * View * Projection).
/// Params:
/// 	fov = The field of view, in radians.
/// 	aspectRatio = The aspect ratio for the viewport this matrix is being used on.
/// 	nearPlane = The closest distance that anything will be rendered.
/// 	viewDistance = The maximum distance to view, such that farPlane is equal to nearPlane + viewDistance.		
Matrix4f fieldOfView(float fov, float aspectRatio, float nearPlane, float viewDistance) {
	float farPlane = nearPlane + viewDistance;
	assert(farPlane > nearPlane && nearPlane > 0 && fov > 0);								
	Matrix4f result = void;
	float calcFov = cast(float)(1f / tan(fov * 0.5f));
	float fovAR = calcFov / aspectRatio;
	result.m11 = fovAR;
	result.m12 = result.m13 = result.m14 = 0;
	result.m22 = calcFov;
	result.m21 = result.m23 = result.m24 = 0;
	result.m31 = result.m32 = 0;
	result.m33 = farPlane / (nearPlane - farPlane);
	result.m34 = -1;
	result.m41 = result.m42 = result.m44 = 0;
	result.m43 = (nearPlane * farPlane) / (nearPlane - farPlane);
	return result;		
}

/// creates a matrix used to look at the given target from the given position.
/// This Matrix is valid as both a DirectX and OpenGL View Matrix, but for OpenGL must be transposed.
/// Note that glUniformMatrix4fv contains a transpose parameter.
/// Params:
/// 	position = The position of the eye, in world space.
/// 	target = The position of the target, in world space.
/// 	up = A vector representing the up direction from the position.
Matrix4f lookAt(Vector3f position, Vector3f target, Vector3f up) {		
	auto first = (position - target).normalized;
	auto second = up.cross(first).normalized;
	auto third = first.cross(second);
	auto result = Matrix4f.identity;
	result.m11 = second.x;
	result.m12 = third.x;
	result.m13 = first.x;
	result.m21 = second.y;
	result.m22 = third.y;
	result.m23 = first.y;
	result.m31 = second.z;
	result.m32 = third.z;
	result.m33 = first.z;
	result.m41 = -second.dot(position);
	result.m42 = -third.dot(position);
	result.m43 = -first.dot(position);
	result.m44 = 1;			
	return result;
}

/// Creates a Matrix to apply a rotation around the X, Y, or Z axis.
/// Params:
/// 	radians = The amount of radians to rotate around the axis.
Matrix4f rotateX(float radians) {
	Matrix4f result = Matrix4f.identity;
	float X = cast(float)cos(radians);
	float Y = cast(float)sin(radians);
	result.m22 = X;
	result.m23 = Y;
	result.m32 = -Y;
	result.m33 = X;
	return result;
}

/// Ditto
Matrix4f rotateY(float radians) {
	Matrix4f result = Matrix4f.identity;
	float X = cast(float)cos(radians);
	float Y = cast(float)sin(radians);
	result.m11 = X;
	result.m13 = -Y;
	result.m31 = Y;
	result.m33 = X;
	return result;
}

/// Ditto
Matrix4f rotateZ(float radians) {
	Matrix4f result = Matrix4f.identity;
	float X = cast(float)cos(radians);
	float Y = cast(float)sin(radians);
	result.m11 = X;
	result.m12 = Y;
	result.m21 = -Y;
	result.m22 = X;
	return result;
}

/// Creates a scalar matrix with the given scale.
Matrix4f scale(float scale) {
	Matrix4f result = Matrix4f.identity;
	result.m11 = scale;
	result.m22 = scale;
	result.m33 = scale;		
	result.m44 = 1;
	return result;	
}

/// Ditto
Matrix4f scale(in Vector3f scale...) {
	Matrix4f result = Matrix4f.identity;
	result.m11 = scale.x;
	result.m22 = scale.y;
	result.m33 = scale.z;
	result.m44 = 1;
	return result;
}

/// Creates a Matrix with the given translation value.
Matrix4f translate(in Vector3f translation ...) {
	Matrix4f result = Matrix4f.identity;
	result.translation = translation;
	return result;
}

@name("Camera Matrix Tests")
unittest {
	auto m1 = fieldOfView(PI_2, 4/3f, 0.1f, 100f);
	auto m2 = lookAt(Vector3f(0, 0, 1), Vector3f(1, 1, 1), Vector3f(0, 1, 0));
	auto rotated = rotateX(PI_2) * rotateY(PI_4) * rotateZ(PI_2);
	auto scaled = scale(2);
	assert(scaled == scale(Vector3f(2, 2, 2)));
	auto translated = translate(Vector3f(1, 1, 1));
	assert(translated.translation == Vector3f(1, 1, 1));
}

alias Matrix!(2, float) Matrix2f;
alias Matrix!(2, double) Matrix2d;
alias Matrix!(2, int) Matrix2i;
alias Matrix!(3, float) Matrix3f;
alias Matrix!(3, double) Matrix3d;
alias Matrix!(3, int) Matrix3i;
alias Matrix!(4, float) Matrix4f;
alias Matrix!(4, double) Matrix4d;
alias Matrix!(4, int) Matrix4i;