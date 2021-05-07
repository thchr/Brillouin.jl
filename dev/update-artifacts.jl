# Deployment of parsed SeeK's data to a Github release as a tar ball, for easy
# integration with Julia's Artifact system

using Pkg.Artifacts, LibGit2, ghr_jll

if !haskey(ENV, "GITHUB_TOKEN")
    error("For automatic github deployment, export GITHUB_TOKEN!")
end

data_dir = @__DIR__

# Where we will put our tarballs
tempdir = mktempdir()

function get_git_remote_url(repo_path::String)
    repo = LibGit2.GitRepo(repo_path)
    origin = LibGit2.get(LibGit2.GitRemote, repo, "origin")
    return LibGit2.url(origin)
end

# Try to detect where we should upload these weights to (or just override
# as shown in the commented-out line)
origin_url = get_git_remote_url(dirname(@__DIR__))
deploy_repo = "$(basename(dirname(origin_url)))/$(splitext(basename(origin_url))[1])"

tag = "data_v0.1.0"
name = "data-SeeK"

@info("Generating $name artifact")
# Create a local artifact
hash = create_artifact() do artifact_dir
    # Copy in parsed data file (TODO: create via `parse-SeeK.jl`)
    cp(joinpath(@__DIR__, "..", "assets", "data", "data-SeeK.jl"), 
       joinpath(artifact_dir, "data-SeeK.jl"))
end

# Spit tarballs to be hosted out to local temporary directory:
tarball_hash = archive_artifact(hash, joinpath(tempdir, "$(name).tar.gz"))

# Calculate tarball url
tarball_url = "https://github.com/$(deploy_repo)/releases/download/$(tag)/$(name).tar.gz"

# Bind this to an Artifacts.toml file
@info("Binding $name in Artifacts.toml...")
bind_artifact!(joinpath(@__DIR__, "..", "Artifacts.toml"), name, hash; download_info=[(tarball_url, tarball_hash)], lazy=true, force=true)

# Upload tarballs to a special github release
@info("Uploading tarballs to $(deploy_repo) tag `$(tag)`")
ghr() do ghr_exe
    run(`$ghr_exe -t -replace -u $(dirname(deploy_repo)) -r $(basename(deploy_repo)) $(tag) $(tempdir)`)
end

@info("Artifacts.toml file now contains all bound artifact names")
