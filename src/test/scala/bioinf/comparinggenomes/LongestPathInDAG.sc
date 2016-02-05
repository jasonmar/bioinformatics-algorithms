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
import bioinf.ComparingGenomes._
import bioinf.Answer._
import bioinf.Input._
import scala.io.Source
val it = Source.fromFile("""C:\tmp\dataset_245_7.txt""").getLines()
val source = it.next().toInt
val sink = it.next().toInt
val edges = readDAG(it.toIndexedSeq)
val answer = longestPathInDAG(source,sink,edges)
//write answer to disk
import java.nio.file.{Paths, Files}
import java.nio.charset.StandardCharsets
Files.write(Paths.get("c:\\tmp\\results.txt"),answer.getBytes(StandardCharsets.UTF_8))
