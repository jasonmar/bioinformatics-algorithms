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

object HiddenMessages {

  /**
    * CODE CHALLENGE: Implement PatternCount (reproduced below).
    * Input: Strings Text and Pattern.
    * Output: Count(Text, Pattern).

    * PatternCount(Text, Pattern)
    * count ← 0
    * for i ← 0 to |Text| − |Pattern|
    * if Text(i, |Pattern|) = Pattern
    * count ← count + 1
    * return count

    * Sample Input:
    * GCGCG
    * GCG

    * Sample Output:
    * 2
    */
  def patternCount(text:String, pattern:String): Int = {
    var count = 0
    for (i <- 0 to text.length - pattern.length) {
      if (text.substring(i,i + pattern.length) == pattern)
        count += 1
    }
    count
  }

  /**
    * CODE CHALLENGE: Solve the Frequent Words Problem.
    * Input: A string Text and an integer k.
    * Output: All most frequent k-mers in Text.

    * Sample Input:
    * ACGTTGCATGTCGCATGATGCATGAGAGCT
    * 4

    * Sample Output:
    * CATG GCAT
    */
  def frequentWords(text:String, k:Int): IndexedSeq[String] = {
    val kmerCounts = kmerFrequency(text,k)
    val maxCount = kmerCounts.values.max // n for most common kmer(s)
    kmerCounts.filter(_._2 == maxCount).keys.toIndexedSeq.sorted // collection of most common kmers
  }

  // extracted for use in other algorithms
  // This function examines text for substrings of length k
  // it returns a map which allows lookup of the number of times a given kmer was found in the text
  @inline
  private def kmerFrequency(text:String,k:Int): Map[String,Int] = {
    getKmers(text,k).groupBy(s => s).map(t => (t._1,t._2.length)) // Map of (kmer,n)
  }

  @inline
  def getKmers(text:String,k:Int): IndexedSeq[String] = {
    val kmers = new Array[String](text.length - k + 1) // preallocate for performance instead of using map
    for (i <- 0 to text.length - k) {
      kmers(i) = text.substring(i, i + k)
    }
    kmers.toIndexedSeq
  }

  @inline
  private def complement(p:Char): Char = {
    p.toUpper match {
      case 'A' => 'T'
      case 'C' => 'G'
      case 'G' => 'C'
      case 'T' => 'A'
      case _ => 'N'
    }
  }

  /**
    Given a nucleotide p, we denote its complementary nucleotide as p. The reverse complement of a string Pattern = p1…pn is the string Pattern = pn … p1 formed by taking the complement of each nucleotide in Pattern, then reversing the resulting string. We will need the solution to the following problem throughout this chapter:

    Reverse Complement Problem: Find the reverse complement of a DNA string.
         Input: A DNA string Pattern.
         Output: Pattern, the reverse complement of Pattern.

    CODE CHALLENGE: Solve the Reverse Complement Problem.

    Sample Input:
         AAAACCCGGT

    Sample Output:
         ACCGGGTTTT

    */
  def reverseComplement(pattern: String): String = {
    val nucleotides = CharBuffer.allocate(pattern.length)
    for (i <- pattern.indices.reverse) {
      nucleotides.put(complement(pattern.charAt(i)))
    }
    new String(nucleotides.array())
  }

  /**
    CODE CHALLENGE: Solve the Pattern Matching Problem.
         Input: Two strings, Pattern and Genome.
         Output: A collection of space-separated integers specifying all starting positions where Pattern appears
         as a substring of Genome.

    Sample Input:
         ATAT
         GATATATGCATATACTT

    Sample Output:
         1 3 9
    */
  def patternMatch(pattern:String,genome:String): IndexedSeq[Int] = {
    (0 to genome.length - pattern.length).flatMap{i =>
      genome.substring(i,i + pattern.length) match {
        case s if s == pattern => Some(i)
        case _ => None
      }
    }
  }

  /**
    CODE CHALLENGE: Solve the Clump Finding Problem (restated below).
    You will need to make sure that your algorithm is efficient enough to handle a large dataset.

    Clump Finding Problem: Find patterns forming clumps in a string.
         Input: A string Genome, and integers k, L, and t.
           k size of k-mer
           L size of window
           t frequency threshold
         Output: All distinct k-mers forming (L, t)-clumps in Genome.

    You can solve the Clump Finding Problem by applying your algorithm for the
    Frequent Words Problem to each window of length L in Genome.

    Sample Input:
         CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
         5 50 4

    Sample Output:
         CGACA GAAGA
    */
  def clumpFindingNaive(genome:String,k:Int,L:Int,t:Int): IndexedSeq[String] = {
    // this implementation does not use dynamic programming and is expected to be slow for large input
    // 150 second runtime for |genome| = 9128 k = 9 L = 598 t = 19
    (0 to genome.length - L).flatMap{i => // combine results for each window
      kmerFrequency(genome,k) // get
        .filter(_._2 >= t).keys
    }.sorted.distinct
  }
}
