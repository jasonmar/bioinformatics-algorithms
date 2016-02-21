/*
 * Bioinformatics Algorithms in Scala
 * Copyright (C) 2016  Jason Mar
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bioinf

import scala.StringBuilder
import scala.collection.mutable
import scala.math.max
import bioinf.Input.finalizeMatrix

object ComparingGenomes {

  /**
    * CODE CHALLENGE: Solve the Change Problem. The DPCHANGE pseudocode is reproduced below for your convenience.
    * Input: An integer money and an array Coins = (coin1, ..., coind).
    * Output: The minimum number of coins with denominations Coins that changes money.

    * Sample Input:
    * 40
    * 50,25,20,10,5,1

    * Sample Output:
    * 2

    * DPCHANGE(money, Coins)
    * MinNumCoins(0) ← 0
    * for m ← 1 to money
    * MinNumCoins(m) ← ∞
    * for i ← 1 to |Coins|
    * if m ≥ coini
    * if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
    * MinNumCoins(m) ← MinNumCoins(m - coini) + 1
    * output MinNumCoins(money)

    */
  def dpChange(money: Int, coins: IndexedSeq[Int]): Int = {
    val minNumCoins = getMinNumCoins(money,coins)
    minNumCoins(money)
  }

  def getMinNumCoins(money:Int, coins: IndexedSeq[Int]): IndexedSeq[Int] = {
    val minNumCoins = new Array[Int](money + 1)
    minNumCoins(0) = 0
    for (m <- 1 to money) {
      minNumCoins(m) = Int.MaxValue
      for (i <- coins.indices) {
        val coin = coins(i)
        if (m >= coin) {
          if (minNumCoins(m - coin) + 1 < minNumCoins(m)) {
            minNumCoins(m) = minNumCoins(m - coin) + 1
          }
        }
      }
    }
    minNumCoins.toIndexedSeq
  }

  /**
    * CODE CHALLENGE: Find the length of a longest path in the Manhattan Tourist Problem.
    * Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right.
    * The two matrices are separated by the - symbol.
    * Output: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid
    * whose edges are defined by the matrices Down and Right.

    * MANHATTANTOURIST(n, m, Down, Right)
    * s0, 0 ← 0
    * for i ← 1 to n
    * si, 0 ← si-1, 0 + downi - 1, 0
    * for j ← 1 to m
    * s0, j ← s0, j−1 + right0, j - 1
    * for i ← 1 to n
    * for j ← 1 to m
    * si, j ← max{si - 1, j + downi - 1, j, si, j - 1 + righti, j - 1}
    * return sn, m

    * Sample Input:
    * 4 4
    * 1 0 2 4 3  Down   n x m+1
    * 4 6 5 2 1
    * 4 4 5 2 1
    * 5 6 8 5 3
    * -
    * 3 2 4 0    Right  n+1 x m
    * 3 2 4 2
    * 0 7 3 3
    * 3 3 0 2
    * 1 3 2 2

    * Sample Output:
    * 34
    */
  def longestPath(n: Int, m: Int, down: IndexedSeq[IndexedSeq[Int]], right: IndexedSeq[IndexedSeq[Int]]): Int = {
    val s = Array.fill[Int](n+1,m+1) { 0 } // initialize to zero

    // s(y axis)(x axis)
    // s(n)(m)
    for (i <- 1 to n) {// solve movement down by simple cumulative sum of downward movement weights
      s(i)(0) = s(i-1)(0) + down(i - 1)(0)
    }
    for (j <- 1 to m) {// solve movement right by simple cumulative sum of rightward movement weights
      s(0)(j) = s(0)(j-1) + right(0)(j - 1)
    }

    for (i <- 1 to n) {
      for (j <- 1 to m) {
        s(i)(j) = max(s(i-1)(j) + down(i-1)(j), s(i)(j-1) + right(i)(j-1))
      }
    }

    s(n)(m)
  }

  /**
    * CODE CHALLENGE: Use OUTPUTLCS (reproduced below) to solve the Longest Common Subsequence Problem.
    * Input: Two strings s and t.
    * Output: A longest common subsequence of s and t. (Note: more than one solution may exist,
    * in which case you may output any one.)
    *
    * OUTPUTLCS(backtrack, v, i, j)
    * if i = 0 or j = 0
    *   return
    * if backtracki, j = "↓"
    *   OUTPUTLCS(backtrack, v, i - 1, j)
    * else if backtracki, j = "→"
    *   OUTPUTLCS(backtrack, v, i, j - 1)
    * else
    *   OUTPUTLCS(backtrack, v, i - 1, j - 1)
    *   output vi
    *
    * Sample Input:
    *   AACCTTGG
    *   ACACTGTGA

    * Sample Output:
    *   AACTGG
    */
  def longestCommonSubsequence(s: String,t: String): String = {
    val n = s.length
    val m = t.length
    val backtrack = LCSBacktrack(s, t)
    val sb = new StringBuilder(n + m)
    outputLCS(backtrack, s, n - 1, m - 1, sb)
    sb.mkString
  }

  def alignmentGraph(v: String, w: String): IndexedSeq[IndexedSeq[Int]] = {
    val n = v.length
    val m = w.length
    val s = Array.fill[Int](n,m){0}
    for (i <- 0 until n) {
      for (j <- 0 until m) {
        // if there is a match, set the score at the current location to 1
        if (v.charAt(i) == w.charAt(j)) s(i)(j) = 1
      }
    }
    finalizeMatrix(s)
  }

  val DIAG = 1 // 1 - match/mismatch \
  val RIGHT = 2 // 2 - insertion -
  val DOWN = 3 // 3 - deletion |

  def LCSBacktrack(v: String, w: String): IndexedSeq[IndexedSeq[Int]] = {
    val n = v.length
    val m = w.length
    val s = LCSPaths(v,w)
    val backtrack = Array.fill[Int](n,m){0}

    if (s(0)(0) == 1) backtrack(0)(0) = DIAG
    for (i <- 1 until n) {backtrack(i)(0) = DOWN}
    for (j <- 1 until m) {backtrack(0)(j) = RIGHT}

    for (i <- 1 until n) {
      for (j <- 1 until m) {
        backtrack(i)(j) = s(i)(j) match {
          case x if x == s(i-1)(j) => DOWN
          case x if x == s(i)(j-1) => RIGHT
          case _ => DIAG
        }
      }
    }
    finalizeMatrix(backtrack)
  }

  def LCSPaths(v: String, w: String): IndexedSeq[IndexedSeq[Int]] = {
    import scala.math.max
    val n = v.length
    val m = w.length
    val s = Array.fill[Int](n,m){0}

    if (v.charAt(0) == w.charAt(0)) s(0)(0) = 1
    for (i <- 1 until n) {s(i)(0) = s(i-1)(0)}
    for (j <- 1 until m) {s(0)(j) = s(0)(j-1)}

    for (i <- 1 until n) {
      for (j <- 1 until m) {
        val charMatch: Boolean = v.charAt(i) == w.charAt(j)
        val nodeAbove = s(i-1)(j)
        val nodeLeft = s(i)(j-1)
        val nodeDiag = {
          if (charMatch) s(i - 1)(j - 1) + 1
          else 0
        }
        // Assign value to current node
        s(i)(j) = max(nodeDiag,max(nodeAbove,nodeLeft))
      }
    }
    finalizeMatrix(s)
  }

  /**
    * STOP and Think: OUTPUTLCS is a recursive algorithm, but it is efficient.
    * What makes it different from the inefficient recursive algorithms for
    * making change and finding a longest path in a DAG?
    * Answer: It always moves towards (0,0) so recursion depth is finite
    *
    * Scala specific: a StringBuilder is used to store the output
    * An implicit could be used to make the invocation look exactly as printed
    */
  def outputLCS(backtrack: IndexedSeq[IndexedSeq[Int]], v: String, i: Int, j: Int, sb: StringBuilder): Unit = {
    if (i >= 0 && j >= 0) {
      backtrack(i)(j) match {
        case DOWN => outputLCS(backtrack, v, i - 1, j, sb) // move UP
        case RIGHT => outputLCS(backtrack, v, i, j - 1, sb) // move LEFT
        case DIAG => outputLCS(backtrack, v, i - 1, j - 1, sb) // move DIAGONAL
          sb.append(v.charAt(i))
        case _ =>
      }
    }
  }

  /**
    * CODE CHALLENGE: Solve the Longest Path in a DAG Problem.
    * Input: An integer representing the source node of a graph, followed by an integer representing the
    * sink node of the graph, followed by a list of edges in the graph. The edge notation 0->1:7 indicates
    * that an edge connects node 0 to node 1 with weight 7.
    * Output: The length of a longest path in the graph, followed by a longest path.
    * (If multiple longest paths exist, you may return any one.)
    *
    * Sample Input:
    * 0
    * 4
    * 0->1:7
    * 0->2:4
    * 2->3:2
    * 1->4:1
    * 3->4:3
    *
    * Sample Output:
    * 9
    * 0->2->3->4
    */
  def longestPathInDAG(source: Int, sink: Int, edges: IndexedSeq[IndexedSeq[Int]]): String = {
    val s = edges.filter(_(0) >= source).filter(_(1) <= sink).sortBy(_(0)) // remove irrelevant edges
    val distance = Array.fill[Int](sink + 1){Int.MinValue}
    val prev = Array.fill[Int](sink + 1){Int.MinValue}
    distance(source) = 0
    val outboundEdges = s.indices.groupBy(i => s(i)(0)) // lookup by source

    def follow(i: Int): Unit = {
      val node0 = s(i)(0)
      val node1 = s(i)(1)
      val weight = s(i)(2)
      val distance1 = distance(node0) + weight
      if (distance(node1) < distance1) {
        distance(node1) = distance1 // keep greatest distance
        prev(node1) = node0
      }
    }

    (source to sink).foreach{node =>
      outboundEdges.getOrElse(node, IndexedSeq[Int]()).foreach(follow) // get all outbound edges for a node, then follow them
    }

    val path = mutable.Stack[Int]()
    path.push(sink)
    var previousNode = sink

    while (previousNode > source){
      previousNode = prev(previousNode)
      path.push(previousNode)
    }

    val sb = new StringBuilder(128)
    sb.append(distance(sink).toString)
    sb.append("\n")
    while (path.nonEmpty) {
      sb.append(path.pop())
      if (path.nonEmpty) sb.append("->")
    }
    sb.mkString
  }


  /**
      *CODE CHALLENGE: Solve the Global Alignment Problem.
         *Input: Two protein strings written in the single-letter amino acid alphabet.
         *Output: The maximum alignment score of these strings followed by an alignment achieving this
         *maximum score. Use the BLOSUM62 scoring matrix and indel penalty σ = 5.

      *Download BLOSUM62 scoring matrix

      *Sample Input:
           *PLEASANTLY
           *MEANLY

      *Sample Output:
           *8
           *PLEASANTLY
           *-MEA--N-LY
    */
  def globalAlignment(v: String, w: String, sigma: Int): GlobalAlignmentResult = {
    val n = v.length + 1
    val m = w.length + 1
    val s = Array.fill[Int](n, m){0} // Alignment Graph

    s(0)(0) = 0 // origin
    for (i <- 1 until n){ // left column
      val j = 0
      s(i)(j) = s(i-1)(j) - sigma
    }
    for (j <- 1 until m){ // top row
      val i = 0
      s(i)(j) = s(i)(j - 1) - sigma
    }

    for (i <- 1 until n) {
      for (j <- 1 until m) {
        val mu = BLOSUM62score(v.charAt(i-1), w.charAt(j-1))
        s(i)(j) = IndexedSeq[Int](
          s(i-1)(j) - sigma, // deletion penalty
          s(i)(j-1) - sigma, // insertion penalty
          s(i-1)(j-1) + mu // match or mismatch value from BLOSUM62 scoring matrix
        ).max
      }
    }
    val alignmentMatrix = finalizeMatrix(s)
    GlobalAlignmentResult(v, w, alignmentMatrix, sigma)
  }

  def scoreGlobalAlignment(x: GlobalAlignmentResult): Int = {
    val s = x.print.split(System.lineSeparator())
    val v = s(1)
    val w = s(2)
    scoreGlobalAlignment(v,w)
  }

  def scoreGlobalAlignment(v:String, w:String): Int = {
    var score = 0
    v.indices.foreach{i =>
      if (v.charAt(i) == '-' || w.charAt(i) == '-') score -= 5
      else {
        score += BLOSUM62score(v.charAt(i), w.charAt(i))
      }
    }
    score
  }

  /** Example Alignment Matrix:
    * 8
    * PLEASANTLY
    * -MEA--N-LY
    * 		     M	 E	 A	 N	 L	 Y
    *      0	 -5	-10	-15	-20	-25	-30
    *  P	-5	 -2	 -6	-11	-16	-21	-26
    *  L	-10	 -3	 -5	 -7	-12	-12	-17
    *  E	-15	 -8	  2	 -3	 -7	-12	-14
    *  A	-20	-13	 -3	  6	  1	 -4	 -9
    *  S	-25	-18	 -8	  1	  7	  2	 -3
    *  A	-30	-23	-13	 -4	  2	  6	  1
    *  N	-35	-28	-18	 -9	  2	  1	  4
    *  T	-40	-33	-23	-14	 -3	  1	 -1
    *  L	-45	-38	-28	-19	 -8	  1	  0
    *  Y	-50	-43	-33	-24	-13	 -4	  8
    *
    * 		  M	  E	   A	  N	  L	  Y
    *  P   (0) -5	 -10	-15	-20	-25	-30  -
    *  L	(-5) -2	  -6	-11	-16	-21	-26  M
    *  E	-10	(-3)  -5	 -7	-12	-12	-17  E
    *  A	-15	 -8	  (2)  -3	 -7	-12	-14  A
    *  S	-20	-13	  -3	 (6)  1	 -4	 -9  -
    *  A	-25	-18	  -8	 (1)	7	  2	 -3  -
    *  N	-30	-23	 -13	(-4)  2	  6	  1  N
    *  T	-35	-28	 -18	 -9	 (2)	1	  4  -
    *  L	-40	-33	 -23	-14	(-3)  1	 -1  L
    *  Y	-45	-38	 -28	-19	 -8	 (1)	0  Y
    *   	-50	-43	 -33	-24	-13	 -4	 (8)
    *
    */
  case class GlobalAlignmentResult(v: String, w: String, alignmentMatrix: IndexedSeq[IndexedSeq[Int]], sigma: Int) {
    def n = v.length + 1
    def m = w.length + 1
    def score = alignmentMatrix(n-1)(m-1) // find the location with greatest alignment score
    lazy val backtrackMatrix = {
      val a = alignmentMatrix
      val s = Array.fill[Int](n, m){0} // backtrackMatrix
      for (i <- 1 until n){s(i)(0) = DOWN} // column 0
      for (j <- 1 until m){s(0)(j) = RIGHT} // row 0
      for (i <- 1 until n){
        for (j <- 1 until m){
          if (a(i)(j) == a(i-1)(j-1) + BLOSUM62score(v.charAt(i-1),w.charAt(j-1))){ // match or mismatch
            s(i)(j) = DIAG
          } else if (a(i)(j) == a(i-1)(j) - sigma){ // deletion
            s(i)(j) = DOWN
          } else if (a(i)(j) == a(i)(j-1) - sigma){ // insertion
            s(i)(j) = RIGHT
          }
        }
      }
      finalizeMatrix(s)
    }
    def print = {
      val len = n + m
      val sbv = new StringBuilder(len)
      val sbw = new StringBuilder(len)

      var i = n-1
      var j = m-1
      var count = 0
      val s = backtrackMatrix
      while (i > 0 || j > 0 || count > len) {
        val path = s(i)(j)
        assert(path == DIAG || path == DOWN || path == RIGHT)
        if (path == DIAG){ // match or mismatch
            assert(i >= 0)
            assert(j >= 0)
            sbv.append(v.charAt(i-1))
            sbw.append(w.charAt(j-1))
            i -= 1
            j -= 1
        } else if (path == DOWN){ // deletion
            assert(i >= 0)
            sbv.append(v.charAt(i-1))
            sbw.append("-")
            i -= 1
        } else if (path == RIGHT){ // insertion
            assert(j >= 0)
            sbv.append("-")
            sbw.append(w.charAt(j-1))
            j -= 1
        }
        count += 1
      }
      val vAlignment = sbv.reverseContents().mkString
      val wAlignment = sbw.reverseContents().mkString
      val result = new StringBuilder(len * 2)
      result.append(score)
      result.append(System.lineSeparator())
      result.append(vAlignment)
      result.append(System.lineSeparator())
      result.append(wAlignment)
      result.mkString
    }
  }

  def BLOSUM62score(v: Char, w: Char): Int = {
    val iv = bs62id(v)
    val iw = bs62id(w)
    BLOSUM62(iv)(iw)
  }

  val bs62id = {
    val chars = IndexedSeq('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
    chars.indices.map{i => (chars(i),i)}.toMap
  }

  val BLOSUM62 = IndexedSeq[IndexedSeq[Int]](
    IndexedSeq(4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2),
    IndexedSeq(0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2),
    IndexedSeq(-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3),
    IndexedSeq(-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2),
    IndexedSeq(-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3),
    IndexedSeq(0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3),
    IndexedSeq(-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2),
    IndexedSeq(-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1),
    IndexedSeq(-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2),
    IndexedSeq(-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1),
    IndexedSeq(-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1),
    IndexedSeq(-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2),
    IndexedSeq(-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3),
    IndexedSeq(-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1),
    IndexedSeq(-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2),
    IndexedSeq(1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2),
    IndexedSeq(0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2),
    IndexedSeq(0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1),
    IndexedSeq(-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2),
    IndexedSeq(-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7)
  )



}
