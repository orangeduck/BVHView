$exeDir = (Get-Item -Path ".\").Parent.FullName
$exeName = "bvhview.exe"
$exePath = $exeDir + "\" + $exeName

$argBvhPath = "../assets/samples/trn_2023_v0_057_main-agent.bvh"
$argWavPath = "../assets/samples/trn_2023_v0_057_main-agent.wav"
$arguments = "-f `"$argBvhPath`" -p `"$argWavPath`""

Start-Process -FilePath $exePath -ArgumentList $arguments
