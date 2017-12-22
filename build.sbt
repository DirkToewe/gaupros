import org.scalajs.jsenv.nodejs.NodeJSEnv
import sbtcrossproject.{CrossType, crossProject}

// enablePlugins(ScalaJSPlugin)

name := "GauProS root project"
scalaVersion in ThisBuild := "2.12.4"

// scalaJSStage in Global := FullOptStage
scalaJSStage in Global := FastOptStage

scalacOptions ++= Seq("-feature", "-deprecation")

lazy val root = project.in( file(".") )
  .aggregate(gauprosJS, gauprosJVM)
  .settings()

lazy val gaupros = crossProject(JSPlatform, JVMPlatform)
  .crossType(CrossType.Pure)
  .settings(
    name := "GauProS",
    version := "0.2",
    scalaVersion := "2.12.4",

    testFrameworks += new TestFramework("utest.runner.Framework"),

    libraryDependencies += "com.lihaoyi" %%% "utest" % "0.4.8" % "test"
  )
  .jsSettings(
    jsEnv := new NodeJSEnv( NodeJSEnv.Config().withArgs("--max_old_space_size=8192" :: Nil) ),
    scalaJSUseMainModuleInitializer := true
  )
  .jvmSettings(

  )
  .in( file(".") )

lazy val gauprosJS = gaupros.js
lazy val gauprosJVM= gaupros.jvm