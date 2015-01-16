package abeellab.macaw

import net.sf.samtools.SAMRecordIterator
import net.sf.samtools.BAMFileReader
import scala.collection.JavaConversions._
import net.sf.samtools.SAMFileReader
import java.io.File
import be.abeel.util.CountMap
import be.abeel.util.TimeInterval
import java.io.PrintWriter
import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult
import net.sf.jannot.utils.SequenceTools
import net.sf.jannot.refseq.Sequence
import net.sf.jannot.refseq.MemorySequence
import scala.io.Source
import atk.util.Tool

import java.text.NumberFormat
import java.util.Locale

import be.abeel.util.FrequencyMap
import be.abeel.util.FrequencyMapUtils

/**
 *
 * SNP typer
 * Features:
 * - 100% matches only
 * -runs entire BAM file
 *
 */
object MacawSNPtyper extends Tool {

  override val version="""
    2015/01/16:    Initial release
    
    """
  
  def revcomp(read: Array[Byte]) = {
    val out = Array.ofDim[Byte](read.length)
    for (i <- 0 until read.length) {
      out(read.length - i - 1) = SequenceTools.complement(read(i).toChar).toByte
    }
    out
  }
  case class Config(val markerFile: String = null, val outputFile: String = null, files: List[File] = List(),val threshold:Int=5)
  /**
   * args(0) = output file
   *
   *
   */
  def main(args: Array[String]): Unit = {

   
    val parser = new scopt.OptionParser[Config]("java -jar macaw.jar") {
      opt[String]("marker") action { (x, c) => c.copy(markerFile = x) } text ("File containing marker sequences. This file has to be a multi-fasta file with the headers indicating the name of the markers.") //, { v: String => config.spacerFile = v })
      opt[String]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]('t',"threshold")action { (x, c) => c.copy(threshold = x) } text ("Threshold to determine absence or presence of a marker (default=5)")
      
      arg[File]("<file>...") unbounded () required () action { (x, c) => c.copy(files = c.files :+ x) } text ("input files")

    }
    parser.parse(args, Config()) map { config =>
      /* Load spacers */
      val lines = if (config.markerFile != null) tLines(config.markerFile).toList  else scala.io.Source.fromInputStream(MacawSNPtyper.getClass().getResourceAsStream("/subset3LongMarkers.txt")).getLines().toList;
      val pw = if (config.outputFile != null) new PrintWriter(config.outputFile) else new PrintWriter(System.out)

      pw.println(generatorInfo)

      pw.println("# Marker input: " + config.markerFile)

      val in = (lines.grouped(2).map(f => (f(0).substring(1), f(1))).toList)

      val repeatSequences = List(("repeat", "GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGA"), ("left_repeat", "GTTTCCGTCCCC"), ("middle_repeat", "TCTCGGGGTTTT"), ("right_repeat", "GGGTCTGACGA"))

      val forwardSpacers = in ++ repeatSequences // ++ mirus()

      val rcSpacers = forwardSpacers.map { case (x, y) => ("RC_" + x, new String(revcomp(y.getBytes))) }

      val spacers = forwardSpacers ++ rcSpacers

    
      pw.println("# Input files: " + config.files)

      /* Prep index */
      val tree = new AhoCorasick();
      for ((id, spacer) <- spacers) {
        tree.add(spacer.getBytes(), id)
      }
      tree.prepare();

      val nf = NumberFormat.getInstance(Locale.US)
      nf.setMaximumFractionDigits(6)

      for (inputFile <- config.files) {
        pw.println("# Processing: " + inputFile)
        /* Connect to bam file*/
        val inputSam = new SAMFileReader(inputFile);

        val cm = new CountMap[String]

        /* Iterate over bamfile */
        val it: SAMRecordIterator = inputSam.iterator()
        var progress: Int = 0
        var totalCoverage: Long = 0

        val time = System.currentTimeMillis();
        val nano = System.nanoTime()

        while (it.hasNext()) {
          val sr = it.next()

          val read = sr.getReadBases()

          totalCoverage += sr.getReadLength()

          val result = tree.search(read);

          for (a <- result) {
            for (s <- a.asInstanceOf[SearchResult].getOutputs()) {
              cm.count(s.asInstanceOf[String])
            }
          }

          progress += 1
          if (progress % 100000 == 0) {
            print(".")


          }
        }
        println

        /* close bam file */
        it.close
        val listx = spacers.filter(p => !p._1.contains("repeat")).map(f => cm.get(f._1).toInt)

        pw.println("# number of spacers = " + spacers.size)
        
        pw.println("## == SNP typing ==")
        for ((x, y) <- spacers) {
          val z = cm.get(x)
          pw.println("#" + x + "\t" + y + "\t" + z)
        }

        val groupedSpacers = spacers.groupBy(pair => pair._1.replaceAll("RC_", ""))

        println("GS: " + groupedSpacers.mkString("\n"))
        val buffer = new StringBuffer()
        pw.println("# Marker\tread-depth\tp-value\tA/P")
        var idx = 0

        println("KS: " + cm.keySet())
        for (gs <- groupedSpacers.filterNot(_._1.contains("repeat")).toList.sortBy(_._1)) {
          idx += 1
          assert(gs._2.size == 2)
          println("GGS: " + gs)
          val z1 = (gs._2.map(p => cm.get(p._1).toInt)).toList
          val z = z1.sum


          pw.println(gs._1 + "\t" + z + "\t" + nf.format(if (z >= config.threshold) 0 else 1) + "\t" + (if (z >= config.threshold) "present" else "absent"))
          buffer.append(if (z >= config.threshold) "1" else "0")
          if (idx % 10 == 0)
            buffer.append(" ")

        }
        pw.println("## Digital markertype: \n" + buffer.toString().grouped(10).mkString(" "))
        pw.println()

      }
      pw.close
    } getOrElse {
      println("Could not interpret command-line arguments, quitting!")
      System.exit(-1)
    }

  }
}
