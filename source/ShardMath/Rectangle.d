module ShardMath.Rectangle;
import ShardMath.Vector;
import std.traits;
import tested;

@nogc:

/// Represents a rectangle containing an x coordinate, y coordinate, width, and height.
struct Rectangle(T) {

	alias Vector!(2, T) Point;
	
	/// Instantiates a new instance of the Rectangle structure.
	/// Params:
	///		x = The x coordinate for this Rectangle.
	///		y = The y coordinate for this Rectangle.
	///		width = The width of this Rectangle.
	///		height = The height of this Rectangle.
	this(T x, T y, T width, T height) {
		this.x = x; 
		this.y = y;
		this.width = width; 
		this.height = height;
	}
	
	/// Returns the smallest Rectangle capable of containing the specified amount of points.
	/// Params: points = The array of points to create a Rectangle from.
	static Rectangle!(T) fromPoints(in Point[] points ...) {
		if(points.length == 0)
			return Rectangle();
		T minX = T.max, minY = T.max;
		T maxX = minVal!T, maxY = minVal!T;
		foreach(ref point; points) {
			minX = point.x < minX ? point.x : minX;
			minY = point.y < minY ? point.y : minY;
			maxX = point.x > maxX ? point.x : maxX;
			maxY = point.y > maxY ? point.y : maxY;
		}		
		return Rectangle(minX, minY, maxX - minX, maxY - minY);
	}

	private static T minVal(T)() pure nothrow @safe {
		static if(isFloatingPoint!T)
			return -T.max;
		else
			return T.min;
	}

	/// Returns an empty rectangle.
	@property static Rectangle!(T) empty() {
		return Rectangle!(T)(0, 0, 0, 0);
	}
	
	/*/// Determines how this Rectangle contains the specified Rectangle.
	/// BUG: This method is not fully implemented yet.
	@disable
	ContainmentType Contains(const ref Rectangle!(T) other) const {
		if(x >= other.x || y >= other.y)
			return ContainmentType.Disjoints;
		if(right >= other.right && bottom >= other.bottom)
			return ContainmentType.Contains;
		return ContainmentType.Intersects;		
	}*/
	
	/// Determines whether this Rectangle contains the specified point.
	bool contains(in Point point) const {
		return x <= point.x && point.x <= (x + width) && y <= point.y && point.y <= (y + height);
	}
	
	/// Returns a Point containing the position of this Rectangle.
	@property Point position() const {
		return Point(x, y);
	}	
	
	/// Returns a Point containing the size of this Rectangle.
	@property Point size() const {
		return Point(width, height);
	}
	
	/// Returns the right-most coordinate in this rectangle.
	@property T right() const {
		return x + width;
	}
	
	/// Returns the bottom-most coordinate in this Rectangle.
	@property T bottom() const {
		return y + height;		
	}

	const Type opCast(Type)() if(hasMember!(Type, "x") && hasMember!(Type, "y") && hasMember!(Type, "width") && hasMember!(Type, "height")) {	
		Type result;
		result.x = cast(typeof(result.x))x;
		result.y = cast(typeof(result.y))y;
		result.width = cast(typeof(result.width))width;
		result.height = cast(typeof(result.height))height;
		return result;
	}
	
	union {
		struct {
			/// The top-left coordinate for this rectangle.
			T x;
			/// The top-right coordinate for this rectangle.
			T y;
			/// The width of this rectangle.
			T width;
			/// The height of this rectangle.
			T height;
		}
		/// Provides array access to the elements of this rectangle.
		T[4] elements;
	}
}

@name("Rectangle Tests")
unittest {
	auto r1 = Rectanglef(1, 1, 2, 3);
	assert(r1.elements == cast(float[4])[1, 1, 2, 3]);
	assert((cast(Rectanglei)r1).elements == cast(int[4])[1, 1, 2, 3]);
	assert(r1.position == Vector2f(1, 1));
	assert(r1.size == Vector2f(2, 3));
	assert(r1.right == 3);
	assert(r1.bottom == 4);
	assert(r1.contains(Vector2f(1, 1)));
	assert(r1.contains(Vector2f(2, 3)));
	assert(!r1.contains(Vector2f(4, 4)));
	assert(!r1.contains(Vector2f(2, 5)));
	assert(Rectanglef.fromPoints(Vector2f(1, 1), Vector2f(3, 4)) == r1);
}

alias Rectangle!(int) Rectanglei;
alias Rectangle!(float) Rectanglef; 