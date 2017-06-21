package markerQuant {

trait ActionObject {

  val description = "A description"

  def main(args: Array[String]) = {
    println("A main function")
  }

  def usage = {
    println("Something to help people")
  }

}

}
