package abeellab.macaw

import atk.util.Tool
import net.sf.jannot.utils.SequenceTools
import java.io.File
import net.sf.jannot.source.FileSource
import scala.collection.JavaConversions._
import be.abeel.util.FrequencyMap


object FindTargetMarkers extends Tool {

  override val version = """
    2015/05/22:    Initial version
    """

  def revcomp(read: Array[Byte]) = {
    val out = Array.ofDim[Byte](read.length)
    for (i <- 0 until read.length) {
      out(read.length - i - 1) = SequenceTools.complement(read(i).toChar).toByte
    }
    out
  }
  case class Config(val target: File = null, val global: File = null, val outputFile: File = null, val markerLength: Integer = 21, val count: Integer = 10)
  /**
   * args(0) = output file
   *
   *
   */
  def main(args: Array[String]): Unit = {

    val parser = new scopt.OptionParser[Config]("java -jar macaw.jar") {
      opt[File]("target") required () action { (x, c) => c.copy(target = x) } text ("File containing target sequences to detect")
      opt[File]("global") required () action { (x, c) => c.copy(global = x) } text ("Reference genome to check uniquely")
      opt[File]('o', "output") action { (x, c) => c.copy(outputFile = x) } text ("File where you want the output to be written")
      opt[Int]("length") action { (x, c) => c.copy(markerLength = x) } text ("Length of markers (default=21nt)")
      opt[Int]("count") action { (x, c) => c.copy(count = x) } text ("Number of markers per target (default=10)")

    }
    parser.parse(args, Config()) map { config =>

      val k=config.markerLength
      
      /*  Make hashset of reference */
      //assume fasta format
      val esTarget=new FileSource(config.target).read()
      println(esTarget.toList)
      
      
      
      /*  Make hashset of targets */
      val esGlobal=new FileSource(config.global).read()
      println(esGlobal.toList)
      
      /*
       * Build index of all kmers with their count in the global sequence
       * 
       * Go with window over each target, 
       * create each 1 bp mismatch sequence and 
       * check all for presence in global, 
       * report if sum of all == 1 
       *  */
      val fm=new FrequencyMap
      
      for(targetTY<-esTarget.toList){
        val target=targetTY.sequence().stringRepresentation()
        for(i<-0 until target.size-k){
          val h=target.substring(i,i+k)
        }
      }
      
      
      
      

    }

  }

}