; Inno Setup script for the minimol Windows installer.
; Build after PyInstaller (dist\minimol must exist):
;   iscc /DAppVersion=0.4.0 packaging\windows\minimol.iss
; Installs for the current user only (no administrator rights needed).

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

[Setup]
AppId={{8F7B6B2E-3C1D-4E3A-9B57-5C2E1A0D4F11}
AppName=minimol
AppVersion={#AppVersion}
AppVerName=minimol {#AppVersion}
AppPublisherURL=https://github.com/fizzest/minimol
AppSupportURL=https://github.com/fizzest/minimol/issues
DefaultDirName={autopf}\minimol
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir=..\..\dist
OutputBaseFilename=minimol-{#AppVersion}-setup
SetupIconFile=minimol.ico
UninstallDisplayIcon={app}\minimol.exe
UninstallDisplayName=minimol
LicenseFile=..\..\LICENSE
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
ChangesAssociations=yes

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked
Name: "openwith"; Description: "Offer minimol in ""Open with"" for .xyz, .mol, .sdf, .mol2, .pdb and .cdxml files"; GroupDescription: "File types:"; Flags: unchecked

[Files]
Source: "..\..\dist\minimol\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\minimol"; Filename: "{app}\minimol.exe"
Name: "{autodesktop}\minimol"; Filename: "{app}\minimol.exe"; Tasks: desktopicon

[Registry]
Root: HKA; Subkey: "Software\Classes\minimol.molecule"; ValueType: string; ValueName: ""; ValueData: "Molecule file"; Flags: uninsdeletekey; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\minimol.molecule\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\minimol.exe,0"; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\minimol.molecule\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\minimol.exe"" ""%1"""; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.xyz\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.mol\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.sdf\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.mol2\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.pdb\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith
Root: HKA; Subkey: "Software\Classes\.cdxml\OpenWithProgids"; ValueType: string; ValueName: "minimol.molecule"; ValueData: ""; Flags: uninsdeletevalue; Tasks: openwith

[Run]
Filename: "{app}\minimol.exe"; Description: "{cm:LaunchProgram,minimol}"; Flags: nowait postinstall skipifsilent
