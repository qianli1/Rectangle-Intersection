#include <array>
#include <chrono>
#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <algorithm>
#include <string_view>
#include <map>

namespace rec2d {

struct Interval {
    double low, high;
    int index;
};
    
    
struct Rectangle {
    const double left;
    const double right;
    const double bottom;
    const double top;
    Rectangle(const std::array<double, 2> &x, const std::array<double, 2> &y)
        : left(x[0]), right(x[1]), bottom(y[0]), top(y[1]) {
        if (left > right || bottom > top)
            throw std::runtime_error("Invalid rectangle");
    }
};
    
struct Point {
    double value;
    int isLeft; //-1: left, 1: right;
    int index;
    
    Point (double _value, int _isLeft, int _index)
        :value(_value), isLeft(_isLeft), index(_index) {
            if (isLeft != -1 && isLeft != 1) {
                throw std::runtime_error("Invalid endpoint");
            }
        }
};
  
  bool operator<(const Point &a, const Point &b) {
    return a.value < b.value;
  }
    
  bool checkOverlap (const Interval &i1, const Interval &i2) {
      if (i1.low > i2.high) return false;
      if (i2.low > i1.high) return false;
      return true;
  }
  
  
  class IntervalIntersections {
    std::map<int, Interval> _inserted;
    public:
    IntervalIntersections(){}
    void insert(const Interval& interval) {_inserted[interval.index]=interval;}
    void remove(const Interval& interval) {
      _inserted.erase(interval.index);
    }
    void find(Interval query, std::vector<std::array<int,2>> &out){
        for (auto &mp: _inserted) {
             if (checkOverlap(query, mp.second)) {
                 out.push_back({query.index, mp.first});
             }
        }
    }
  };

//struct IntervalTreeNode {
//    Interval *interval;
//    double max;
//    IntervalTreeNode *left, *right;
//};

class IntervalTreeNode {
private:
    Interval * interval;
    double max;
    IntervalTreeNode *left, *right;
    
    IntervalTreeNode * newNode(const Interval &i);
    IntervalTreeNode * minLowInterval();
    static double findMax(double a,double b,double c);
    
public:
    IntervalTreeNode(){};
    IntervalTreeNode* getLeft() {return left;};
    IntervalTreeNode* insert(const Interval& interval);
    IntervalTreeNode* remove(const Interval& interval);
    void find(Interval query, std::vector<std::array<int,2>> &out);
};

double IntervalTreeNode::findMax(double a, double b, double c){
    return a > b ? (a > c ? a : c) : (b > c ? b : c);
}

IntervalTreeNode * IntervalTreeNode::minLowInterval()
{
    IntervalTreeNode* current = left;
 
    while (current->getLeft() != NULL)
        current = current->getLeft();
 
    return current;
}

IntervalTreeNode* IntervalTreeNode::newNode(const Interval &i){
    IntervalTreeNode *temp = new IntervalTreeNode;
    temp->interval = new Interval(i);
    temp->max = i.high;
    temp->left = temp->right = NULL;
    return temp;
}

IntervalTreeNode* IntervalTreeNode::insert(const Interval &i) {
    if (interval == NULL)
        return newNode(i);

    // Get low value of interval at root
    double l = interval->low;
    
    // If root's low value is smaller, then new interval goes to left subtree
    if (i.low < l)
        left = left->insert(i);

    // Else, new node goes to right subtree.
    else
        right = right->insert(i);

    // Update the max value of this ancestor if needed
    if (max < i.high)
        max = i.high;

    return this;
    
}

IntervalTreeNode * IntervalTreeNode::remove(const Interval& i) {
    if(interval==NULL){
        return NULL;
    }
    if(i.low < interval->low){
        left=left->remove(i);
    }
    else if(i.low > interval->low){
        right=right->remove(i);
    }
    else if(i.low == interval->low){
        if(i.high == interval->high){
            if(left == NULL){
                IntervalTreeNode *temp = right;
                free(this);
                return temp;
            }
            else if(right == NULL){
                IntervalTreeNode *temp = left;
                free(this);
                return temp;
            }
            IntervalTreeNode *temp = right->minLowInterval();
            interval = temp->interval;
            right = right->remove(*(temp->interval));
            }
        else
            right = right->remove(i);
        
    }
    max=findMax(interval->high,((left!=NULL)?left->max:LONG_MIN),((right!=NULL)?right->max:LONG_MIN));
    return this;
}

void IntervalTreeNode::find(Interval query, std::vector<std::array<int,2>> &out) {
    if (interval == NULL) return;

    // If given interval overlaps with root
    if (checkOverlap(*(interval), query))
    {
        std::array<int, 2> overlapPair = {query.index, interval->index};
        out.push_back(overlapPair);
    }
    
    if (left != NULL && left->max >= query.low)
       left->find(query, out);

    right->find(query, out);
}
 
template <class F>
std::pair<std::vector<std::array<int, 2>>, std::chrono::milliseconds> benchmark(F f, const std::vector<Rectangle> &rs) {
    auto start = std::chrono::steady_clock::now();
    auto ct = f(rs);
    auto end = std::chrono::steady_clock::now();
    return {ct, std::chrono::duration_cast<std::chrono::milliseconds>(end - start)};
}

std::string getWKT(const std::vector<Rectangle> &rs) {
    std::ostringstream out;
    out << "MULTIPOLYGON(";
    for (std::size_t i = 0; i < rs.size(); ++i) {
        if (i != 0)
            out << ", ";
        out << "((";
        out << rs[i].left << " " << rs[i].bottom << ", ";
        out << rs[i].right << " " << rs[i].bottom << ", ";
        out << rs[i].right << " " << rs[i].top << ", ";
        out << rs[i].left << " " << rs[i].top << ", ";
        out << rs[i].left << " " << rs[i].bottom << ", ";
        out << "))";
    }
    out << ")";
    return out.str();
}

std::vector<rec2d::Rectangle> getRandomRectangles(const Rectangle &maxBounds, const std::array<double, 2> &maxExtents,
                                                 std::size_t n) {
    auto gen = std::mt19937(std::random_device()());
    std::uniform_real_distribution<double> xd(maxBounds.left, maxBounds.right);
    std::uniform_real_distribution<double> yd(maxBounds.bottom, maxBounds.top);
    std::uniform_real_distribution<double> xed(0, maxExtents[0]);
    std::uniform_real_distribution<double> yed(0, maxExtents[1]);

    std::vector<rec2d::Rectangle> ret;
    for (std::size_t i = 0; i < n; ++i) {
        auto xl = xd(gen);
        auto yl = yd(gen);
        ret.emplace_back(Rectangle{{xl, xl + xed(gen)}, {yl, yl + yed(gen)}});
    }
    return ret;
}

bool intersectPairRectangle (const Rectangle & r1, const Rectangle & r2) {
    if (r1.left > r2.right || r1.right < r2.left) return false;
    if (r1.bottom > r2.top || r1.top < r2.bottom) return false;
    return true;
}
    
// compute the a list of indexes of intersecting pairs of rectangles, in a simple, clear manner
std::vector<std::array<int, 2>> getIntersectionsSimple(const std::vector<Rectangle> &rs) {
    std::vector<std::array<int, 2>> simple;
    for (std::size_t i = 0; i < rs.size(); ++i) {
        for (std::size_t j = i + 1; j < rs.size(); ++j) {
            if (intersectPairRectangle(rs[i], rs[j])) {
                simple.push_back({static_cast<int>(i), static_cast<int>(j)});
            }
        }
    }
    return simple;
}

// let's make it faster
std::vector<std::array<int, 2>> getIntersectionsFast(const std::vector<Rectangle> &rs) {
    std::vector<std::array<int, 2>> fast;
    
    std::vector<Point> points;
    for (size_t i = 0; i < rs.size(); ++i) {
        auto value = rs[i].left;
        //int isLeft = -1;
        points.emplace_back(Point({value, -1, static_cast<int>(i)}));
        value = rs[i].right;
        //isLeft = 1;
        points.emplace_back(Point({value, 1, static_cast<int>(i)}));
    }
    
    std::sort(points.begin(), points.end());
    
    IntervalIntersections tree;
    
    for (auto point : points) {
        auto index = point.index;
        Interval interval = Interval{rs[index].bottom, rs[index].top, index};
        if (point.isLeft == -1) {
            tree.find(interval, fast);
            tree.insert(interval);
        } else if (point.isLeft == 1) {
            tree.remove(interval);
        }
    }
  
         return fast;
}

std::vector<std::array<int, 2>> getIntersectionsFaster(const std::vector<Rectangle> &rs) {
    std::vector<std::array<int, 2>> faster;
    
    std::vector<Point> points;
    for (size_t i = 0; i < rs.size(); ++i) {
        auto value = rs[i].left;
        //int isLeft = -1;
        points.emplace_back(Point({value, -1, static_cast<int>(i)}));
        value = rs[i].right;
        //isLeft = 1;
        points.emplace_back(Point({value, 1, static_cast<int>(i)}));
    }
    
    std::sort(points.begin(), points.end());
    
    IntervalTreeNode tree;
    
    for (auto point : points) {
        auto index = point.index;
        Interval interval = Interval{rs[index].bottom, rs[index].top, index};
        if (point.isLeft == -1) {
            tree.find(interval, faster);
            tree.insert(interval);
        } else if (point.isLeft == 1) {
            tree.remove(interval);
        }
    }
  
         return faster;
}

} // namespace rec2d

int main(int, char **) {
    auto rs = rec2d::getRandomRectangles({{-100., 100.}, {-100., 100.}}, {.1, .1}, 20000);

    auto [faster, timeFaster] = rec2d::benchmark(rec2d::getIntersectionsFast, rs);
    std::cout << "Got " << faster.size() << " in " << timeFaster.count() << " for faster" << std::endl;
    
    auto [fast, timeFast] = rec2d::benchmark(rec2d::getIntersectionsFast, rs);
    std::cout << "Got " << fast.size() << " in " << timeFast.count() << " for fast" << std::endl;

    auto [simple, timeSimple] = rec2d::benchmark(rec2d::getIntersectionsSimple, rs);
    std::cout << "Got " << simple.size() << " in " << timeSimple.count() << " for simple" << std::endl;

    return EXIT_SUCCESS;
}
