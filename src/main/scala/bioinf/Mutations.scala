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

import java.nio.CharBuffer
import bioinf.Input._
import scala.collection.mutable

object Mutations {
  /**
    * CODE CHALLENGE: Solve the Trie Construction Problem.
    * Input: A collection of strings Patterns.
    * Output: The adjacency list corresponding to Trie(Patterns), in the following format.
    * If Trie(Patterns) has n nodes, first label the root with 0 and then label the remaining nodes with
    * the integers 1 through n - 1 in any order you like. Each edge of the adjacency list of
    * Trie(Patterns) will be encoded by a triple: the first two members of the triple must be the
    * integers labeling the initial and terminal nodes of the edge, respectively; the third member
    * of the triple must be the symbol labeling the edge.

    * Sample Input:
    * ATAGA
    * ATC
    * GAT
    * Sample Output:
    * 0->1:A
    * 1->2:T
    * 2->3:A
    * 3->4:G
    * 4->5:A
    * 2->6:C
    * 0->7:G
    * 7->8:A
    * 8->9:T
    */
  def trieConstruction(patterns: IndexedSeq[String]): Trie = {
    val len = patterns.length * patterns.map{_.length}.max

    val node0 = Array.fill[Int](len){-1}
    val node1 = Array.fill[Int](len){-1}
    val label = Array.fill[Int](len){-1}
    val edges = mutable.Map[(Int,Int),Int]() // edge pointers for (nodeId,label) tuples
    var i = 0 // edge pointer
    var j = 0 // node pointer

    patterns.foreach{pattern =>
      var currentNode = 0 // begin at ROOT node
      pattern.map{label2int}.foreach{currentSymbol =>
        val edge = (currentNode,currentSymbol)
        edges.get(edge) match {
          case Some(idx) => // edge exists
            currentNode = node1(idx)
          case _ => // new edge
            j += 1 // increment node pointer
            node0(i) = currentNode // write initial node to adjacency list
            node1(i) = j // write terminal node to adjacency list
            label(i) = currentSymbol // write label to adjacency list
            currentNode = j // continue from node j for next symbol
            edges += ((edge,i)) // add edge to index
            i += 1 // increment edge pointer
        }
      }

    }
    val data = IndexedSeq(node0,node1,label).map(_.slice(0,i).toIndexedSeq)
    Trie(AdjacencyList(data(0),data(1),data(2)))
  }

  case class AdjacencyList(v: IndexedSeq[Int], w: IndexedSeq[Int], label: IndexedSeq[Int])


  /**
    * CODE CHALLENGE: Implement TRIEMATCHING to solve the Multiple Pattern Matching Problem.
    * Input: A string Text and a collection of strings Patterns.
    * Output: All starting positions in Text where a string from Patterns appears as a substring.

    * TRIEMATCHING(Text, Trie)
    *   while Text is nonempty
    *     PREFIXTRIEMATCHING(Text, Trie)
    *     remove first symbol from Text

    * PREFIXTRIEMATCHING(Text, Trie)
    *   symbol ← first letter of Text
    *   v ← root of Trie
    *   while forever
    *     if v is a leaf in Trie
    *       return the pattern spelled by the path from the root to v
    *     else if there is an edge (v, w) in Trie labeled by symbol
    *       symbol ← next letter of Text
    *       v ← w
    *     else
    *       output "no matches found"
    *       return

    * Sample Input:
    *   AATCGGGTTCAATCGGGGT
    *   ATCG
    *   GGGT
    * Sample Output:
    *   1 4 11 15
    */
  def trieMatching(text: String, patterns: IndexedSeq[String]): IndexedSeq[Int] = {
    val buf = CharBuffer.wrap(text).asReadOnlyBuffer()
    val trie = trieConstruction2(patterns)
    val matchIndices = mutable.ArrayBuffer[Int]()

    var symbol = 0

    while (buf.hasRemaining) {
      symbol = label2int(buf.get())
      buf.mark()
      var v = 0
      var continue = true
      while(continue) {
        val edge = trie.edges0.get((v,symbol))
        val leaf = trie.edges0.get((v,-1)) // pattern termination
        if (leaf.nonEmpty) {
          // match is complete
          buf.reset()
          matchIndices += buf.position() - 1
          continue = false
        } else if (edge.nonEmpty && buf.hasRemaining) {
          symbol = label2int(buf.get()) // next letter of text
          v = trie.adjacencyList.w(edge.get) // next node in trie
        } else {
          // no matches found
          buf.reset()
          continue = false
        }
      }
    }
    matchIndices.result().toIndexedSeq
  }

  case class Trie(adjacencyList: AdjacencyList) {
    lazy val leaves = adjacencyList.w.toSet.diff(adjacencyList.v.toSet)
    lazy val edges = adjacencyList.label.indices.map{i => ((adjacencyList.v(i),adjacencyList.w(i)),i)}.toMap
    lazy val edges0 = adjacencyList.label.indices.map{i => ((adjacencyList.v(i),adjacencyList.label(i)),i)}.toMap
    lazy val edges1 = adjacencyList.label.indices.map{i => ((adjacencyList.w(i),adjacencyList.label(i)),i)}.toMap
  }

  /** This trie construction method appends '$' to each pattern to indicate
    * This allows successful resolution of patterns which are nested within each other
    * Example in the book is 'pan' and 'pantry'
    * With terminator characters, they become 'pan$' and 'pantry$' such that 'pan' can be found in the trie
    */
  def trieConstruction2(patterns: IndexedSeq[String]): Trie = {
    trieConstruction(patterns.map{_ + "$"})
  }

  /**
      *CODE CHALLENGE: Solve the Suffix Tree Construction Problem.
         *Input: A string Text.
         *Output: The edge labels of SuffixTree(Text). You may return these strings in any order.

      *Sample Input:
        *ATAAATG$
      *Sample Output:
        *AAATG$
        *G$
        *T
        *ATG$
        *TG$
        *A
        *A
        *AAATG$
        *G$
        *T
        *G$
        *$

    *MODIFIEDSUFFIXTRIECONSTRUCTION(Text)
        *Trie ← a graph consisting of a single node root
        *for i ← 0 to |Text| - 1
            *currentNode ← root
            *for j ← i to |Text| - 1
                *currentSymbol ← j-th symbol of Text
                *if there is an outgoing edge from currentNode labeled by currentSymbol
                    *currentNode ← ending node of this edge
                *else
                    *add a new node newNode to Trie
                    *add an edge newEdge connecting currentNode to newNode in Trie
                    *Symbol(newEdge) ← currentSymbol
                    *Position(newEdge) ← j
                    *currentNode ← newNode
            *if currentNode is a leaf in Trie
                *assign label i to this leaf
        *return Trie


    *
    */
  def suffixTrieConstruction(text:String): SuffixTrie = {
    val v = mutable.ArrayBuffer[Int]()
    val w = mutable.ArrayBuffer[Int]()
    val pos = mutable.ArrayBuffer[Int]()
    val label = mutable.ArrayBuffer[Int]()
    var edgeId = 0 // edge pointer
    var nodeId = 0 // node pointer
    val edges = mutable.Map[(Int,Int),Int]() // edge pointers for (nodeId,label) tuples
    val positionLabel = mutable.Map[Int,Int]() // index pointers for edgeId -> position
    var currentNode = 0
    var symbol = 0
    var edge = (0,0)
    val TERMINATOR_SYMBOL = -1

    for (i <- 0 until text.length) {
      currentNode = 0
      for (j <- i until text.length) {
        symbol = label2int(text.charAt(j))
        edge = (currentNode,symbol)
        edges.get(edge) match {
          case Some(idx) =>
            currentNode = w(idx)
          case _ =>
            // add a new node
            v += currentNode
            nodeId += 1
            w += nodeId

            // add an edge
            edges += ((edge,edgeId))
            pos += j
            label += symbol
            edgeId += 1

            currentNode = nodeId
        }

      }
      if (symbol == TERMINATOR_SYMBOL) {
        positionLabel += ((edges.getOrElse(edge,0),i))
      }
    }
    val data = IndexedSeq(v,w,label,pos).map{_.result().toIndexedSeq}
    SuffixTrie(AdjacencyList(data(0),data(1),data(2)),data(3),positionLabel.toMap)
  }
  case class SuffixTrie(adjacencyList: AdjacencyList, position: IndexedSeq[Int], label: Map[Int,Int]) {
    def indices = adjacencyList.v.indices
    def length = adjacencyList.v.length
    lazy val outEdges = indices.groupBy(i => adjacencyList.v(i))
  }
  case class Path(v: Int, w: Int, pos: Int, len: Int)

  /**
      *MODIFIEDSUFFIXTREECONSTRUCTION(Text)
          *Trie ← MODIFIEDSUFFIXTRIECONSTRUCTION
          *for each non-branching path Path in Trie
              *substitute Path by a single edge e connecting the first and last nodes of Path
              *Position(e) ← Position(first edge of Path)
              *Length(e) ← number of edges of Path
          *return Trie
    */
  def suffixTreeConstruction(text:String): SuffixTree = {
    val text1 = if (text.last == '$') text else text + "$"
    val t = suffixTrieConstruction(text1)

    val stack = mutable.Stack[Int]() // used for DFS
    val v = mutable.ArrayBuffer[Int]()
    val w = mutable.ArrayBuffer[Int]()
    val pos = mutable.ArrayBuffer[Int]()
    val len = mutable.ArrayBuffer[Int]()

    val out = t.outEdges
    val currentNode = 0
    out.get(currentNode).foreach(_.foreach(stack.push))

    while (stack.nonEmpty) {
      val i = stack.pop()
      val startNode = t.adjacencyList.v(i)
      val position = t.position(i)

      var length = 1
      var currentNode = t.adjacencyList.w(i)
      var next: Option[Int] = t.outEdges.get(currentNode) match {
        case Some(ids) if ids.length == 1 => // Path
          Some(ids.head)
        case Some(ids) if ids.length > 1 => // Branch
          ids.foreach(stack.push)
          None
        case _ => // leaf
          None
      }
      while (next.nonEmpty) {
        currentNode = t.adjacencyList.w(next.get)
        length += 1
        next = t.outEdges.get(currentNode) match {
          case Some(ids) if ids.length == 1 => // Path
            Some(ids.head)
          case Some(ids) if ids.length > 1 => // Branch
            ids.foreach(stack.push)
            None
          case _ => // leaf
            None
        }
      }
      v += startNode
      w += currentNode
      pos += position
      len += length
    }

    SuffixTree(text1,Edges(v.toIndexedSeq, w.toIndexedSeq, pos.toIndexedSeq, len.toIndexedSeq))
  }
  case class Edges(v: IndexedSeq[Int], w: IndexedSeq[Int], pos: IndexedSeq[Int], len: IndexedSeq[Int])

  case class SuffixTree(text:String, edges: Edges) {
    lazy val out = edges.v.indices.groupBy(i => edges.v(i))
    lazy val in = edges.v.indices.map{i => (edges.w(i),i)}.toMap
  }

  // returns indices of edges terminating at deepest nodes
  def deepestEdges(s: SuffixTree): IndexedSeq[Int] = {
    val depth = Array.fill[Int](s.edges.v.length){-1}
    val stack = mutable.Stack[Int]() // used for DFS

    var currentNode = 0
    var parentNode = -1
    var maxDepth = -1
    val deepEdgeIds = mutable.Set[Int]()

    s.out.get(currentNode).foreach(_.foreach(stack.push))
    while (stack.nonEmpty) {
      val i = stack.pop()
      parentNode = s.edges.v(i)
      currentNode = s.edges.w(i)
      val currentDepth = s.in.get(parentNode) match {
        case Some(idx) => depth(idx) + s.edges.len(i)
        case _ => s.edges.len(i)
      }
      depth(i) = currentDepth
      s.out.get(currentNode) match {
        case Some(indices) =>
          if (indices.length > 1) { // only return edges that branch
            if (currentDepth > maxDepth) {
              maxDepth = currentDepth
              deepEdgeIds.clear()
              deepEdgeIds += i
            } else if (currentDepth == maxDepth) {
              deepEdgeIds += i
            }
          }
          indices.foreach(stack.push)
        case _ =>
      }
    }
    deepEdgeIds.toIndexedSeq
  }

  // prints the string terminating at edge with given index
  def printSuffix(s: SuffixTree, startIdx: Int): String = {
    var parentNode = s.edges.v(startIdx)
    val sb = new StringBuilder(65536)
    val stack = mutable.Stack[Int]()
    stack.push(startIdx)
    while (stack.nonEmpty) {
      val i = stack.pop()
      parentNode = s.edges.v(i)
      s.text
        .substring(s.edges.pos(i),s.edges.pos(i) + s.edges.len(i))
        .reverseIterator
        .foreach(sb.append)
      s.in.get(parentNode).foreach(stack.push)
    }
    sb.reverseContents().mkString
  }
}
